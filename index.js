
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
        this.m_dens.fill(0);
        this.m_densPrev.fill(0);
    },
    m_pixels: new Array(g_nGridSize).fill(0)
}

// fluid.m_dens = fluid.m_dens.map(x => Math.random())

// fluid.m_densPrev = fluid.m_densPrev.map(x => Math.random())
console.log(fluid.m_dens)

// helper functions for fluids
const seek = (iWidth, iHeight) => iWidth + (g_nWidth + 2) * iHeight

console.log(g_nWidth, g_nHeight)
// set half of prev to one density, the other to another
for (let i = 1; i <= g_nHeight; ++i) {
    for (let j = 1; j <= g_nWidth; ++j) {
        if (j >= 2.5 *g_nWidth / 6 && j <= 3.5 *g_nWidth / 6 && i >= 2 * g_nHeight/ 6 && i <= 4 * g_nHeight/ 6) {
            // fluid.m_densPrev[seek(j, i)] = 0
            fluid.m_dens[seek(j, i)] = 0
            // fluid.m_u[seek(j,i)] = 0
            // fluid.m_uPrev[seek(j,i)] = 0
        } else {

            // console.log(j, g_nWidth / 2)
            // fluid.m_densPrev[seek(j, i)] = 1
            fluid.m_dens[seek(j, i)] = 1
            // fluid.m_u[seek(j,i)] = 1
            // fluid.m_uPrev[seek(j,i)] = 1
        }
    }
}
for (let i = 1; i <= g_nHeight; ++i) {
    for (let j = 1; j <= g_nWidth; ++j) {
        if (j < g_nWidth / 2) {
            fluid.m_u[seek(j,i)] = 1
            fluid.m_uPrev[seek(j,i)] = 1
        } else {
            fluid.m_u[seek(j,i)] = -1
            fluid.m_uPrev[seek(j,i)] = -1
        }
    }
}

for (let i = 0; i <= g_nHeight+1; ++i) {
    for (let j = 0; j <= g_nWidth+1; ++j) {
        fluid.m_pixels[seek(j, i)] = (i / g_nHeight)
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

            // console.log(x, j)

            let y = Math.min(Math.max(i - dt0_h * v[seek(j, i)], 0.5), N + 0.5)
            // console.log(x, y, j, i, u[seek(j, i)])
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
console.log(g_canvas.width, g_canvas.height)

let stop = true
let iFrameCount = 0
let fps = 30
let startTime, lastTime;

function animate(curTime) {
    if (!startTime) {
        startTime = curTime 
        lastTime = curTime
    }

    const dt = curTime - lastTime
    if(dt < 1000/64 && !stop) return requestAnimationFrame(animate);


    g_ctx.clearRect(0, 0, g_nWidth, g_nHeight)
    const visc = 0
    const diff = 0.0001
    vel_step(g_nWidth, g_nHeight, fluid.m_u, fluid.m_v, fluid.m_uPrev, fluid.m_vPrev, visc, 0.024)
    dens_step(g_nWidth, g_nHeight, fluid.m_dens, fluid.m_densPrev, fluid.m_u, fluid.m_v, diff, 0.024)
    // we need to advect the fluid according to velocity fields
    let temp = new Array(g_nGridSize).fill(0)
    advect(g_nWidth, g_nHeight, 0, temp, fluid.m_pixels, fluid.m_u, fluid.m_v, dt)
    fluid.m_pixels = temp;


    // console.log(dt)
    // console.log(fluid.m_dens)
    // diffuse(g_nWidth, g_nHeight, 0, fluid.m_dens, fluid.m_densPrev, 2, dt)
    // fluid.m_densPrev = fluid.m_dens.map(x => x)


    // console.log(fluid.m_dens)

    // draw densities
    // fluid.m_dens = fluid.m_dens.map(x => Math.random())
    for(let i = 1; i < g_nHeight + 1; ++i) {
        for (let j = 1; j < g_nWidth + 1; ++j) {
            // get color
            let dens = fluid.m_dens[seek(j, i)]
            let vel = fluid.m_u[seek(j, i)]
            let pixel = fluid.m_pixels[seek(j, i)]
            g_ctx.fillStyle = `hsl(${dens / 13}turn, 100%, 50%)`
            // g_ctx.fillStyle = `hsl(${pixel * 0.4 + 0.65}turn, 80%, 65%)`
            // const outof255 = Math.floor(Math.sqrt((vel + 1) / 2) * 255)
            // const outof255 = Math.floor(dens * 255)
            // g_ctx.fillStyle = `rgb(${outof255}, ${outof255}, ${outof255})`
            g_ctx.fillRect(j - 1, i - 1, 1, 1)
        }
    }
    iFrameCount += 1
    lastTime = curTime
    if (!stop && iFrameCount< 2000) {
        // console.log(iFrameCount)
        requestAnimationFrame(animate)
    }else {
        console.log("STOPPED", curTime - startTime)
    }
}

function drawNoStep() {
    g_ctx.clearRect(0, 0, g_nWidth, g_nHeight)
    // draw densities
    // fluid.m_dens = fluid.m_dens.map(x => Math.random())
    for(let i = 1; i < g_nHeight + 1; ++i) {
        for (let j = 1; j < g_nWidth + 1; ++j) {
            // get color
            let dens = fluid.m_dens[seek(j, i)]
            let vel = fluid.m_uPrev[seek(j, i)]
            g_ctx.fillStyle = `hsl(${dens / 13}turn, 100%, 50%)`
            // const outof255 = Math.floor(Math.sqrt(vel) * 255)
            // const outof255 = Math.floor(dens * 255)
            // g_ctx.fillStyle = `rgb(${outof255}, ${outof255}, ${outof255})`
            g_ctx.fillRect(j - 1, i - 1, 1, 1)
        }
    }
}

requestAnimationFrame(animate)
// create a timer loop

 
function clickedbutton(e) {
    console.log("hi", e)
    if (stop) {
        e.target.textContent = "Stop"
    } else {
        e.target.textContent = "Start"
    }
    stop = !stop;
    if (!stop) {
        console.log("ANIMATING")
        iFrameCount = 0
        // fluid.m_densPrev = fluid.m_densPrev.map(x => Math.random())
        // fluid.m_dens = fluid.m_densPrev.map(x => Math.random())
        requestAnimationFrame(animate)
    }
}

function resetDensity() {
    for (let i = 0; i <= g_nHeight+1; ++i) {
        for (let j = 0; j <= g_nWidth+1; ++j) {
            if (j >= 2.5 *g_nWidth / 6 && j <= 3.5 *g_nWidth / 6 && i >= 2 * g_nHeight/ 6 && i <= 4 * g_nHeight/ 6) {
                fluid.m_densPrev[seek(j,i)] = 0
            } else {
                fluid.m_densPrev[seek(j,i)] = 1
            }
        }
    }
}

function resetVelocity() {
    for (let i = 0; i <= g_nHeight+1; ++i) {
        for (let j = 0; j <= g_nWidth+1; ++j) {
            if (j >= 2.5 *g_nWidth / 6 && j <= 3.5 *g_nWidth / 6 && i >= 2 * g_nHeight/ 6 && i <= 4 * g_nHeight/ 6) {
                fluid.m_uPrev[seek(j,i)] = 1
                fluid.m_u[seek(j,i)] = 1
            } else {
                fluid.m_uPrev[seek(j,i)] = 0
                fluid.m_u[seek(j,i)] = 0
            }
        }
    }
    console.log("Velocity reset!")
}

function resetAll() {
    fluid.reset()

    // resetDensity()
    resetVelocity()

    iFrameCount = 0;

    requestAnimationFrame(animate)
}

g_canvas.addEventListener("mousedown", (evt) => {
    // get initial coordinates
    const rect = evt.target.getBoundingClientRect();
    const scaleX = g_canvas.width / rect.width;
    const scaleY = g_canvas.height / rect.height;
    let prevX = Math.max(Math.floor((evt.clientX - rect.left) * scaleX), 0) + 1;
    let prevY = Math.max(Math.floor((evt.clientY - rect.top) * scaleY), 0) + 1;
    let prevT = performance.now();
    // add mouse move listener
    let handleMouseMove = (evt) => {
        const x = Math.max(Math.floor((evt.clientX - rect.left) * scaleX), 0) + 1;
        const y = Math.max(Math.floor((evt.clientY - rect.top) * scaleY), 0) + 1;
        const T = performance.now()
        const elapsed = T - prevT
        // set horizontal vel
        const v_x = (x - prevX) / elapsed
        const v_y = (y - prevY) / elapsed
        console.log(v_x, v_y);
        // fluid.m_uPrev = 
        // set vert vel
    }

    g_canvas.addEventListener("mousemove", handleMouseMove)
})

g_canvas.addEventListener("mousemove", (evt) => {
    // console.log(evt);
    const rect = evt.target.getBoundingClientRect();
    const scaleX = g_canvas.width / rect.width;
    const scaleY = g_canvas.height / rect.height;
    const x = Math.max(Math.floor((evt.clientX - rect.left) * scaleX), 0) + 1;
    const y = Math.max(Math.floor((evt.clientY - rect.top) * scaleY), 0) + 1;
    // console.log(x, y)
    fluid.m_densPrev[seek(x, y)] += 100.001
    // evt.offsetX;
    // evt.offsetY;
    requestAnimationFrame(animate)
    // drawNoStep();
})