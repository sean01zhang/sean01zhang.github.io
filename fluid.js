
// Let's try hungarian notation...
const g_canvas = document.getElementById("fluid")
const g_nWidth = g_canvas.width
const g_nHeight = g_canvas.height
const g_nGridSize = (g_nWidth + 2) * (g_nHeight + 2)

// navier stokes arrays
const fluid = {
    m_u: new Array(g_nGridSize).fill(0),
    m_v: new Array(g_nGridSize).fill(0),
    m_uPrev: new Array(g_nGridSize).fill(0),
    m_vPrev: new Array(g_nGridSize).fill(0),
    m_dens: new Array(g_nGridSize).fill(1),
    m_densPrev: new Array(g_nGridSize).fill(1),
    reset: function() {
        this.m_u.fill(0);
        this.m_v.fill(0);
        this.m_uPrev.fill(0);
        this.m_vPrev.fill(0);
        this.m_dens.fill(1);
        this.m_densPrev.fill(1);
    },
    m_pixels: new Array(g_nGridSize).fill(0)
}

// helper functions for fluids
const seek = (iWidth, iHeight) => iWidth + (g_nWidth + 2) * iHeight

// initial pixel state
for (let i = 0; i <= g_nHeight+1; ++i) {
    for (let j = 0; j <= g_nWidth+1; ++j) {
        fluid.m_pixels[seek(j, i)] = (i / g_nHeight)
        fluid.m_uPrev[seek(j, i)] = 2 * Math.random() - 1;
        fluid.m_vPrev[seek(j, i)] = 2 * Math.random() - 1;
    }
}

// create pressure waves
function add_source(M, N, x, s, dt) {
    const size = (M+2) * (N+2)
    for (let i = 0; i < size; ++i) x[i] += dt * s[i]
}

// diffuse (gauss-seidel relaxation solve matrix)
// x0 = x_ij - a(x_{i-1,j} + x_{i+1,j} + x_{i, j-1} + x_{i, j+1} - 4x_ij)
function diffuse(M, N, b, x, x0,diff, dt) {
    const a = dt * diff * M * N

    for (let k = 0; k < 20; ++k) { // I don't quite understand why 20 iterations was chosen.
        for (let i = 1; i <= N; ++i) {
            for (let j = 1; j <= M; ++j) {
                x[seek(j, i)] = (x0[seek(j, i)] + a * (x[seek(j-1, i)] + x[seek(j+1, i)] + x[seek(j, i-1)] + x[seek(j, i+1)])) / (1 + 4*a)
            }
        }
        set_bnd(M, N, b, x)
    }
}

// density advection to some velocity field
function advect(M, N, b, d, d0, u, v, dt) {
    // console.log(u, v, M, N, d.length, u.length, v.length)
    let dt0_h = dt * N / 2
    let dt0_w = dt * M / 2
    for (let i = 1; i <= N; ++i) {
        for (let j = 1; j <= M; ++j) {
            let x = Math.min(Math.max(j - dt0_w * u[seek(j, i)], 0.5), M + 0.5)
            let y = Math.min(Math.max(i - dt0_h * v[seek(j, i)], 0.5), N + 0.5)

            let j0 = Math.floor(x)
            let j1 = j0 + 1
            let i0 = Math.floor(y)
            let i1 = i0 + 1

            let s1 = x - j0
            let s0 = 1 - s1
            let t1 = y - i0
            let t0 = 1 - t1

            d[seek(j, i)] = s0 * (t0 * d0[seek(j0, i0)] + t1 * d0[seek(j0, i1)])
                                + s1 * (t0 * d0[seek(j1, i0)] + t1 * d0[seek(j1, i1)])
        }
    }
    set_bnd(M, N, b, d)
}

function dens_step(M, N, x, x0, u, v, diff, dt) {
    // adding new stuff from x0 to x
    add_source(M, N, x, x0, dt);
    // x -> m_densPrev; x0 -> m_dens
    [x, x0] = [x0, x]
    // m_densPrev now has diffused data
    diffuse(M, N, 0, x, x0, diff, dt);
    // x -> m_dens; x0 -> m_densPrev
    [x, x0] = [x0, x]
    // for (let i = 0; i < g_nGridSize; ++i) {
    //     x[i] = x0[i]
    // }
    advect(M, N, 0, x, x0, u, v, dt)

    // reset x0
    x0.fill(0);
}
// evolving velocities

function vel_step(M, N, u, v, u0, v0, visc, dt) {
    add_source(M, N, u, u0, dt)
    add_source(M, N, v, v0, dt);
    [u, u0] = [u0, u];
    diffuse(M, N, 1, u, u0, visc, dt);
    [v, v0] = [v0, v];
    diffuse(M, N, 2, v, v0, visc, dt);

    project(M, N, u, v, u0, v0);
    [u, u0] = [u0, u];
    [v, v0] = [v0, v];

    advect(M, N, 1, u, u0, u0, v0, dt);
    advect(M, N, 2, v, v0, u0, v0, dt);

    project(M, N, u, v, u0, v0)
    u0.fill(0);
    v0.fill(0);
}

function project(M, N, u, v, p, div) {
    let h = 1 / N
    for (let i = 1; i <= N; ++i) {
        for (let j = 1; j <= M; ++j) {
            div[seek(j, i)] = -0.5*h*(u[seek(j+1, i)] - u[seek(j-1, i)] + v[seek(j, i+1)] - v[seek(j, i-1)])
            p[seek(j, i)] = 0
        }
    }
    set_bnd(M, N, 0, div);
    set_bnd(M, N, 0, p);

    for(let k=0; k < 20; ++k) {
        for (let i = 1; i <= N; ++i) {
            for (let j = 1; j <= M; ++j) {
                p[seek(j, i)] = 1/4 * (div[seek(j, i)] + p[seek(j-1, i)] + p[seek(j+1, i)] + p[seek(j, i-1)] + p[seek(j, i+1)])
            }
        }
        set_bnd(M, N, 0, p)
    }

    for(let i = 1; i <= N; ++i) {
        for(let j = 1; j <= M; ++j) {
            u[seek(j, i)] -= 0.5 * (p[seek(j+1, i)] - p[seek(j-1, i)]) / h
            v[seek(j, i)] -= 0.5 * (p[seek(j, i+1)] - p[seek(j, i-1)]) / h
        }
    }
    set_bnd(M, N, 1, u);
    set_bnd(M, N, 2, v);
}

function set_bnd(M, N, b, x) {
    for(let i = 1; i <= N; ++i) {
        x[seek(0, i)]   = b == 1 ? -x[seek(1, i)] : x[seek(1, i)]
        x[seek(M+1, i)] = b == 1 ? -x[seek(M, i)] : x[seek(M, i)]
    }
    for (let j = 1; j <= M; ++j) {
        x[seek(j, 0)]   = b == 2 ? -x[seek(j, 1)] : x[seek(j, 1)]
        x[seek(j, N+1)] = b == 2 ? -x[seek(j, N)] : x[seek(j, N)]
    }
    // corners
    x[seek(0, 0)]     = 0.5 * (x[seek(1, 0)]   + x[seek(0, 1)])
    x[seek(0, N+1)]   = 0.5 * (x[seek(1, N+1)] + x[seek(0, N)])
    x[seek(M+1, 0)]   = 0.5 * (x[seek(M, 0)]   + x[seek(M+1, 1)])
    x[seek(M+1, N+1)] = 0.5 * (x[seek(M, N+1)] + x[seek(M+1, N)])
}

const g_ctx = g_canvas.getContext("2d")

let stop = false
let iFrameCount = 0
let fps = 30
let startTime, lastTime;
let g_colorOffset= Math.random();
let g_isDark = false;

if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
    // dark mode
    g_isDark = true;
}
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', event => {
    g_isDark = event.matches;
});

function animate(curTime) {
    if (!startTime) {
        startTime = curTime 
        lastTime = curTime
    }

    const dt = curTime - lastTime
    if(dt < 1000/30 && !stop) return requestAnimationFrame(animate);


    g_ctx.clearRect(0, 0, g_nWidth, g_nHeight)
    const visc = 0.00002
    const diff = 0.12
    // console.log(dt/10);
    vel_step(g_nWidth, g_nHeight, fluid.m_u, fluid.m_v, fluid.m_uPrev, fluid.m_vPrev, visc, dt/100)
    dens_step(g_nWidth, g_nHeight, fluid.m_dens, fluid.m_densPrev, fluid.m_u, fluid.m_v, diff, dt/100)
    // we need to advect the fluid according to velocity fields
    let temp = new Array(g_nGridSize).fill(0)
    advect(g_nWidth, g_nHeight, 0, temp, fluid.m_pixels, fluid.m_u, fluid.m_v, dt/100)
    fluid.m_pixels = temp;

    // draw densities
    // fluid.m_dens = fluid.m_dens.map(x => Math.random())
    for(let i = 1; i < g_nHeight + 1; ++i) {
        for (let j = 1; j < g_nWidth + 1; ++j) {
            // get color
            let pixel = fluid.m_pixels[seek(j, i)]

            let saturation = g_isDark ? 74 : 74;
            let lightness = g_isDark ? 25 : 60;

            g_ctx.fillStyle = `hsl(${pixel * 0.5 + g_colorOffset}turn, ${saturation}%, ${lightness}%)`
            g_ctx.fillRect(j - 1, i - 1, 1, 1)
        }
    }
    iFrameCount += 1
    lastTime = curTime
    if (!stop && iFrameCount< 1000) {
        requestAnimationFrame(animate)
    } else {
        // console.log("STOPPED", curTime - startTime)
    }
}

requestAnimationFrame(animate)

let g_prevX = 0;
let g_prevY = 0;
let g_prevT = performance.now();

let handleMouseMove = (evt) => {
    evt.preventDefault();
    const rect = evt.target.getBoundingClientRect();
    const scaleX = g_canvas.width / rect.width;
    const scaleY = g_canvas.height / rect.height;
    const x = Math.max(Math.floor((evt.clientX - rect.left) * scaleX), 0) + 1;
    const y = Math.max(Math.floor((evt.clientY - rect.top) * scaleY), 0) + 1;
    const T = performance.now()
    const elapsed = T - g_prevT

    if (elapsed < 0.000001) {
        //TODO: figure out why elapsed is sometimes 0
        return;
    }

    // set horizontal vel
    const v_x = (x - g_prevX) / elapsed
    if (v_x) {
        fluid.m_uPrev[seek(x, y)] = Math.max(Math.min(v_x, 2.50), -2.5) * 20;
    }
    // set vert vel
    const v_y = (y - g_prevY) / elapsed
    if (v_y) {
        fluid.m_vPrev[seek(x, y)] = Math.max(Math.min(v_y, 2.5), -2.5) * 20;
    }

    g_prevT = T;
    g_prevX = x;
    g_prevY = y;
}

let handleMouseMoveMobile = (evt) => {
    evt.preventDefault();
    const rect = evt.target.getBoundingClientRect();
    const scaleX = g_canvas.width / rect.width;
    const scaleY = g_canvas.height / rect.height;
    const x = Math.max(Math.floor((evt.touches[0].clientX - rect.left) * scaleX), 0) + 1;
    const y = Math.max(Math.floor((evt.touches[0].clientY - rect.top) * scaleY), 0) + 1;
    const T = performance.now()
    const elapsed = T - g_prevT

    if (elapsed < 0.000001) {
        //TODO: figure out why elapsed is sometimes 0
        return;
    }

    // set horizontal vel
    const v_x = (x - g_prevX) / elapsed
    if (v_x) {
        fluid.m_uPrev[seek(x, y)] = Math.max(Math.min(v_x, 2.50), -2.5) * 20;
    }
    // set vert vel
    const v_y = (y - g_prevY) / elapsed
    if (v_y) {
        fluid.m_vPrev[seek(x, y)] = Math.max(Math.min(v_y, 2.5), -2.5) * 20;
    }

    g_prevT = T;
    g_prevX = x;
    g_prevY = y;
}
 
g_canvas.addEventListener("mouseenter", (evt) => {
    evt.preventDefault();
    // get initial coordinates
    const rect = evt.target.getBoundingClientRect();
    const scaleX = g_canvas.width / rect.width;
    const scaleY = g_canvas.height / rect.height;
    g_prevX = Math.max(Math.floor((evt.clientX - rect.left) * scaleX), 0) + 1;
    g_prevY = Math.max(Math.floor((evt.clientY - rect.top) * scaleY), 0) + 1;
    g_prevT = performance.now();
    // add mouse move listener
    g_canvas.addEventListener("mousemove", handleMouseMove);
})

g_canvas.addEventListener("mouseleave", (evt) => {
    evt.preventDefault();
    g_canvas.removeEventListener("mousemove", handleMouseMove);
})

g_canvas.addEventListener("touchstart", (evt) => {
    evt.preventDefault();
    // get initial coordinates
    const rect = evt.target.getBoundingClientRect();
    const scaleX = g_canvas.width / rect.width;
    const scaleY = g_canvas.height / rect.height;
    console.log(rect, scaleX, scaleY, evt.clientX, evt.clientY, evt)
    g_prevX = Math.max(Math.floor((evt.touches[0].clientX - rect.left) * scaleX), 0) + 1;
    g_prevY = Math.max(Math.floor((evt.touches[0].clientY - rect.top) * scaleY), 0) + 1;
    g_prevT = performance.now();
    console.log(g_prevX, g_prevY, g_prevT)
    // add mouse move listener
    g_canvas.addEventListener("touchmove", handleMouseMoveMobile);
})

g_canvas.addEventListener("touchend", (evt) => {
    evt.preventDefault();
    g_canvas.removeEventListener("touchmove", handleMouseMoveMobile);
})