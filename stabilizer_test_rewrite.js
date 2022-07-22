const stabilizers = [
    // identity
    (points) => {
        return points;
    },

    // interpolation
    (points) => {
        const stabilization = 0.9;
        let floatX = points[0][0];
        let floatY = points[0][1];

        let newPoints = [];
        for(let i = 0; i < points.length; i += 1){
            newPoints.push([floatX, floatY]);

            floatX += (points[i][0] - floatX) * (1 - stabilization);
            floatY += (points[i][1] - floatY) * (1 - stabilization);
        }

        return newPoints;
    },
    
    // bezier
    (points) => {
        const maxDistance = 50;
        const calcPointDist = (p1, p2) => Math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2);
    
        let collectedPoints = [points[0]];
        let prevPoint = points[0]; // point before "currentPoint"
        let currentDist = 0;

        for(let i = 0; i < points.length; i++){
            let point = points[i];
            currentDist += calcPointDist(point, prevPoint);

            if(maxDistance < currentDist) {
                collectedPoints.push(point);
                lastPoint = point;
                lastIndex = i;

                currentDist = 0;
            }

            prevPoint = point;
        }

        let pushedLast = collectedPoints[collectedPoints.length - 1];
        let actualLast = points[points.length - 1];
        if(
            pushedLast[0] != actualLast[0] &&
            pushedLast[1] != actualLast[1]
        ) collectedPoints.push(actualLast);

        let n = collectedPoints.length - 1;

        if(n < 2) return collectedPoints;

        let C = new Array(n).fill(0).map((_, i) => {
            const row = new Array(n).fill(0);
            row[i] = 4;

            if(i < n - 1) row[i + 1] = 1;
            if(i > 0) row[i - 1] = 1;

            return row;
        });

        C[0][0] = 2
        C[n - 1][n - 1] = 7
        C[n - 1][n - 2] = 2

        let P = new Array(n).fill(0).map((_, i) => {
            return [
                2 * (2 * collectedPoints[i][0] + collectedPoints[i + 1][0]), 
                2 * (2 * collectedPoints[i][1] + collectedPoints[i + 1][1]),
            ]
        });

        P[0][0] = collectedPoints[0][0] + 2 * collectedPoints[1][0]
        P[0][1] = collectedPoints[0][1] + 2 * collectedPoints[1][1]
        P[n - 1][0] = 8 * collectedPoints[n - 1][0] + collectedPoints[n][0]
        P[n - 1][1] = 8 * collectedPoints[n - 1][1] + collectedPoints[n][1]

        const A0 = lusolve(C, P.map(pair => pair[0]));
        const A1 = lusolve(C, P.map(pair => pair[1]));

        let B = new Array(n).fill(0);
        B = B.map((_, i) => {
            return [
                2 * collectedPoints[i + 1][0] - A0[i + 1],
                2 * collectedPoints[i + 1][1] - A1[i + 1],
            ];
        });
        B[n - 1] = [
            (A0[n - 1] + collectedPoints[n][0]) / 2,
            (A1[n - 1] + collectedPoints[n][1]) / 2,
        ];

        let newPoints = [];
        for(let i = 0; i < n; i++){
            const leftPoint = collectedPoints[i];
            const curA = [A0[i], A1[i]];
            const curB = B[i];
            const rightPoint = collectedPoints[i + 1];

            const lambda = (a, b, c, d) => (t) => [
                (1-t)**3*a[0] + 3*(1-t)**2*t*b[0] + 3*(1-t)*t**2*c[0] + t**3*d[0],
                (1-t)**3*a[1] + 3*(1-t)**2*t*b[1] + 3*(1-t)*t**2*c[1] + t**3*d[1],
            ];
            const newFunc = lambda(leftPoint, curA, curB, rightPoint);

            if(i != n - 1){
                let curSegPoints = new Array(10).fill(0).map((_, i) => i / (10 - 1));
                newPoints.push(...curSegPoints.map(n => newFunc(n)))
            } else {
                let index = Math.floor(10 * currentDist / maxDistance)
                let curSegPoints = new Array(index).fill(0).map((_, i) => i / (10 - 1));
                newPoints.push(...curSegPoints.map(n => newFunc(n)))
            }
        }

        return newPoints;
    },

    // spring    
    (points) => {
        let prevPoint = points[0];
        
        const k = 0.1, dt = 1, m = 1, z = 0.7;

        let velocity = [0, 0];
        let position = points[0];
        let collectedPoints = [];

        for(let i = 0; i < points.length; i++){
            let currPoint = points[i];

            let dx = [
                currPoint[0] - prevPoint[0],
                currPoint[1] - prevPoint[1]
            ];

            let accel = [
                k * dx[0] / m,
                k * dx[1] / m
            ];
            
            velocity = [
                (velocity[0] + accel[0] * dt) * z,
                (velocity[1] + accel[1] * dt) * z
            ];

            position = [
                position[0] + velocity[0] * dt,
                position[1] + velocity[1] * dt
            ];
            
            collectedPoints.push(position);

            prevPoint = position;
        }

        return collectedPoints;
    },

    // savitzky-golay 5th degree
    (points) => {

        const maxDistance = 3;
        const calcPointDist = (p1, p2) => Math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2);

        let collectedPoints = [points[0]];
        let prevPoint = points[0]; // point before "currentPoint"
        let currentDist = 0;

        for(let i = 0; i < points.length; i++){
            let point = points[i];
            currentDist += calcPointDist(point, prevPoint);

            if(maxDistance < currentDist) {
                collectedPoints.push(point);
                lastPoint = point;
                lastIndex = i;

                currentDist = 0;
            }

            prevPoint = point;
        }

        let pushedLast = collectedPoints[collectedPoints.length - 1];
        let actualLast = points[points.length - 1];
        if(
            pushedLast[0] != actualLast[0] &&
            pushedLast[1] != actualLast[1]
        ) collectedPoints.push(actualLast);
        
        let finalPoints = [];

        for(let i = 2; i < collectedPoints.length - 2; i++){
            let a = collectedPoints[i - 2];
            let b = collectedPoints[i - 1];
            let c = collectedPoints[i];
            let d = collectedPoints[i + 1];
            let e = collectedPoints[i + 2];

            let newPoint = [
                (-3*a[0] + 12*b[0] + 17*c[0] + 12*d[0] + -3*e[0])/35,
                (-3*a[1] + 12*b[1] + 17*c[1] + 12*d[1] + -3*e[1])/35,
            ];

            finalPoints.push(newPoint);
        }

        return finalPoints;
    },

    // savitzky-golay 9th degree
    (points) => {

        const maxDistance = 3;
        const calcPointDist = (p1, p2) => Math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2);
    
        let collectedPoints = [points[0]];
        let prevPoint = points[0]; // point before "currentPoint"
        let currentDist = 0;

        for(let i = 0; i < points.length; i++){
            let point = points[i];
            currentDist += calcPointDist(point, prevPoint);

            if(maxDistance < currentDist) {
                collectedPoints.push(point);
                lastPoint = point;
                lastIndex = i;

                currentDist = 0;
            }

            prevPoint = point;
        }

        let pushedLast = collectedPoints[collectedPoints.length - 1];
        let actualLast = points[points.length - 1];
        if(
            pushedLast[0] != actualLast[0] &&
            pushedLast[1] != actualLast[1]
        ) collectedPoints.push(actualLast);
        
        let finalPoints = [];

        /*for(let i = 2; i < collectedPoints.length - 2; i++){
            let a = collectedPoints[i - 2];
            let b = collectedPoints[i - 1];
            let c = collectedPoints[i];
            let d = collectedPoints[i + 1];
            let e = collectedPoints[i + 2];

            let newPoint = [
                (-3*a[0] + 12*b[0] + 17*c[0] + 12*d[0] + -3*e[0])/35,
                (-3*a[1] + 12*b[1] + 17*c[1] + 12*d[1] + -3*e[1])/35,
            ];

            finalPoints.push(newPoint);
        }*/

        /*for(let i = 3; i < collectedPoints.length - 3; i++){
            let a = collectedPoints[i - 3];
            let b = collectedPoints[i - 2];
            let c = collectedPoints[i - 1];
            let d = collectedPoints[i];
            let e = collectedPoints[i + 1];
            let f = collectedPoints[i + 2];
            let g = collectedPoints[i + 3];

            let newPoint = [
                (-2*a[0] + 3*b[0] + 6*c[0] + 7*d[0] + 6*e[0] + 3*f[0] + -2*g[0])/21,
                (-2*a[1] + 3*b[1] + 6*c[1] + 7*d[1] + 6*e[1] + 3*f[1] + -2*g[1])/21,
            ];

            finalPoints.push(newPoint);
        }*/

        for(let i = 4; i < collectedPoints.length - 4; i++){
            let a = collectedPoints[i - 4];
            let b = collectedPoints[i - 3];
            let c = collectedPoints[i - 2];
            let d = collectedPoints[i - 1];
            let e = collectedPoints[i];
            let f = collectedPoints[i + 1];
            let g = collectedPoints[i + 2];
            let h = collectedPoints[i + 3];
            let j = collectedPoints[i + 4];

            let newPoint = [
                (-21*a[0] + 14*b[0] + 39*c[0] + 54*d[0] + 59*e[0] + 54*f[0] + 39*g[0] + 14*h[0] + -21*j[0])/231,
                (-21*a[1] + 14*b[1] + 39*c[1] + 54*d[1] + 59*e[1] + 54*f[1] + 39*g[1] + 14*h[1] + -21*j[1])/231,
            ];

            finalPoints.push(newPoint);
        }

        return finalPoints;
    }
];

class Engine {
    constructor(canvasContainer, scheme){
        this.canvasContainer = canvasContainer;
        const canvasDimensions = canvasContainer.getBoundingClientRect();
        this.width = canvasDimensions.width | 0;
        this.height = canvasDimensions.height | 0;
        
        this.canvas = document.createElement('canvas');
        this.canvas.width = this.width;
        this.canvas.height = this.height;
        canvasContainer.appendChild(this.canvas);

        this.ctx = this.canvas.getContext('2d');
        this.scheme = scheme;

        this.rawStrokes = [];
        this.finalStrokes = [];
        this.finalPixels = [];
        this.SDF;

        this.penDown = false;
    }

    createStroke(){
        this.rawStrokes.push([]);
        this.finalStrokes.push([]);
        this.finalPixels.push([]);
    }

    updateStroke(point){
        this.rawStrokes[this.rawStrokes.length - 1].push(point);

        let rawPoints = this.rawStrokes[this.rawStrokes.length - 1];
        this.finalStrokes[this.finalStrokes.length - 1] = stabilizers[this.scheme](rawPoints);

        this.clearCanvas();
        this.drawStrokes(referenceStroke);
    }

    drawStrokes(reference = null){
        this.finalStrokes.forEach(stroke => {
            this.ctx.beginPath();

            stroke.forEach((point, i) => {
                if(i == 0) this.ctx.moveTo(point[0], point[1]);
                else this.ctx.lineTo(point[0], point[1]);
            });

            // this.ctx.closePath();
            this.ctx.strokeStyle = 'black';
            this.ctx.stroke();
        });


        if(reference != null){
            this.ctx.beginPath();

            reference.forEach((point, i) => {
                if(i == 0) this.ctx.moveTo(point[0], point[1]);
                else this.ctx.lineTo(point[0], point[1]);
            });

            // this.ctx.closePath();
            this.strokeStyle = 'gray';
            this.ctx.stroke();
        }
        // this.calculateSDF();
    }

    returnPixels(){
        return this.ctx.getImageData(
            0, 0,
            this.width, this.height
        );
    }

    calculateSDF(){
        console.log("Calculating SDF for schema: ", this.scheme)
        const getSDF = () => {
            for(let i = 0; i < this.finalStrokes.length; i++){
                let currentStroke = this.finalStrokes[i];
    
                for(let j = 0; j < currentStroke.length - 1; j++){
                    this.finalPixels[i].push(
                        ...bline(
                            currentStroke[j][0] | 0, 
                            currentStroke[j][1] | 0,
                            currentStroke[j + 1][0] | 0, 
                            currentStroke[j + 1][1] | 0,
                        )
                    );
                }
            }

            let imageData = this.ctx.getImageData(0, 0, this.width, this.height);
            const newImageData = this.ctx.createImageData(imageData);
            const SDF = new Array(this.width * this.height).fill(0);
            const fastDist = (a, b, c, d) => (a - c) ** 2 + (b - d) ** 2;
            let pixels = imageData.data;
            let newPixels = newImageData.data;

            for(let i = 0; i < pixels.length; i += 4){
                let index =  i;
                let r = pixels[index]; 
                let g = pixels[index + 1];
                let b = pixels[index + 2];
                let a = pixels[index + 3];

                if(a != 0){
                    newPixels[index] = r;
                    newPixels[index + 1] = g;
                    newPixels[index + 2] = b;
                    newPixels[index + 3] = 255;
                    
                    SDF[i / 4] = 1;
                } else {
                    let smallestDist = Infinity;
                    
                    let x = (i / 4) % this.width | 0;
                    let y = (i / 4) / this.width | 0;

                    for(let p of this.finalPixels[this.finalPixels.length - 1]){
                        let curDist = fastDist(x, y, p[0], p[1]);
                        if(curDist < smallestDist){
                            smallestDist = curDist;
                        }
                    }

                    let pixelDist = Math.sqrt(smallestDist);
                    let ratio = pixelDist / 100;

                    newPixels[index] = 255 * ratio;
                    newPixels[index + 1] = 255 * ratio;
                    newPixels[index + 2] = 255;
                    newPixels[index + 3] = 255;

                    SDF[i / 4] = ratio;
                }
            }

            this.ctx.putImageData(newImageData, 0, 0);
            this.SDF = SDF;
        }

        getSDF();
    }

    compareSDFV2(pixelReference){
        console.log("Calculating SDF for schema: ", this.scheme)
        
        // calculate SDF
        const calcSDF = (width, height, pixels = []) => {
            const SDF = new Array(width * height).fill(0);
            const fastDist = (a, b, c, d) => (a - c) ** 2 + (b - d) ** 2;

            // for each pixel
            for(let i = 0; i < SDF.length * 4; i += 4){
                let x = (i / 4) % this.width | 0;
                let y = (i / 4) / this.width | 0;

                if(pixels.some(n => n[0] == x & n[1] == y)){
                    SDF[i / 4] = -1;
                } else {
                    // if not on stroke
                    // find distance from stroke
                    let smallestDist = Infinity;
                    
                    for(let p of pixels){
                        let curDist = fastDist(x, y, p[0], p[1]);
                        if(curDist < smallestDist){
                            smallestDist = curDist;
                        }
                    }
    
                    let pixelDist = Math.sqrt(smallestDist);
                    let ratio = pixelDist / 100;
    
                    SDF[i / 4] = ratio;
                }
            }

            return SDF;
        }

        // convert points to pixels
        for(let i = 0; i < this.finalStrokes.length; i++){
            let currentStroke = this.finalStrokes[i];

            for(let j = 0; j < currentStroke.length - 1; j++){
                this.finalPixels[i].push(
                    ...bline(
                        currentStroke[j][0] | 0, 
                        currentStroke[j][1] | 0,
                        currentStroke[j + 1][0] | 0, 
                        currentStroke[j + 1][1] | 0,
                    )
                );
            }
        }

        const newImageData = this.ctx.createImageData(this.width, this.height);
        let newPixels = newImageData.data;
        const SDF = calcSDF(this.width, this.height, this.finalPixels[this.finalPixels.length - 1]);
        const referenceSDF = calcSDF(this.width, this.height, pixelReference);

        // draw reference stroke
        for(let i = 0; i < SDF.length * 4; i += 4){
            if (referenceSDF[i / 4] == -1){
                newPixels[i] = 200;
                newPixels[i + 1] = 200;
                newPixels[i + 2] = 200;
                newPixels[i + 3] = 255;

                continue;
            }
            
            /*if (SDF[i / 4] == -1){
                newPixels[i] = 0;
                newPixels[i + 1] = 0;
                newPixels[i + 2] = 0;
                newPixels[i + 3] = 255;

                continue;
            } else if (referenceSDF[i / 4] == -1){
                newPixels[i] = 255;
                newPixels[i + 1] = 255;
                newPixels[i + 2] = 255;
                newPixels[i + 3] = 255;

                continue;
            }
            newPixels[i] = 255 * (referenceSDF[i / 4]);
            newPixels[i + 1] = 255 * (referenceSDF[i / 4]);
            newPixels[i + 2] = 255;
            newPixels[i + 3] = 255;*/
        }

        // draw user stroke using overlapping SDF data
        console.log(this.finalPixels)
        for(let pixel of this.finalPixels[this.finalPixels.length - 1]){
            let i = (pixel[1] * this.width + pixel[0]) * 4;

            newPixels[i] = 255 * referenceSDF[i / 4] * 8;
            newPixels[i + 1] = 0;
            newPixels[i + 2] = 0;
            newPixels[i + 3] = 255;
        }

        this.ctx.putImageData(newImageData, 0, 0);
        this.SDF = SDF;
    }

    calculateCurvature(){
        const stroke = this.finalStrokes[this.finalStrokes.length - 1];

        let displacements = [];
        let totalLength = 0;
        let prevPt = stroke[0];

        for(let i = 0; i < stroke.length; i++){
            let currentPt = stroke[i];
            let dx = currentPt[0] - prevPt[0];
            let dy = currentPt[1] - prevPt[1];

            displacements[i] = [dx, dy];
            totalLength += Math.sqrt(dx ** 2, dy ** 2);

            prevPt = currentPt;
        }

        let firstDerivative = [];
        let totalFstDDist = 0;
        prevPt = displacements[0];

        for(let i = 0; i < displacements.length; i++){
            let currentPt = displacements[i];
            let dx = currentPt[0] - prevPt[0];
            let dy = currentPt[1] - prevPt[1];

            firstDerivative[i] = [dx, dy];
            totalFstDDist += Math.sqrt(dx ** 2, dy ** 2);

            prevPt = currentPt;
        }

        let secondDerivative = [];
        let totalSndDDist = 0;
        prevPt = firstDerivative[0];

        for(let i = 0; i < firstDerivative.length; i++){
            let currentPt = firstDerivative[i];
            let dx = currentPt[0] - prevPt[0];
            let dy = currentPt[1] - prevPt[1];

            secondDerivative[i] = [dx, dy];
            totalSndDDist += Math.sqrt(dx ** 2, dy ** 2);

            prevPt = currentPt;
        }

        return [totalLength, totalFstDDist, totalSndDDist];
    }

    compareSDF(inputSDF){
        let sum = 0;

        for(let i = 0; i < inputSDF.length; i++){
            sum += Math.abs(this.SDF[i] - inputSDF[i]);
        }

        return sum;
    }

    clearCanvas(){
        this.ctx.clearRect(
            0, 0, 
            this.canvas.width, this.canvas.height
        );
    }
}

const canvasA = document.getElementById('canvasA');
const canvasB = document.getElementById('canvasB');
const canvasC = document.getElementById('canvasC');
const canvasD = document.getElementById('canvasD');
const canvasE = document.getElementById('canvasE');
const canvasF = document.getElementById('canvasF');

//let canvases = [canvasA, canvasB, canvasC, canvasD, canvasE, canvasF];
let canvases = [canvasA,  canvasD, canvasF];
let engines = canvases.map((canvas, i) => new Engine(canvas, i));

let mouseDown = false;
canvasA.addEventListener('pointerdown', () => {
    mouseDown = true;
    engines.map(e => e.penDown = true);
    engines.map(e => e.createStroke());
});
canvasA.addEventListener('pointerup', () => {
    mouseDown = false;
    engines.map(e => e.penDown = false);

    engines.map(e => e.compareSDFV2(referenceStroke));
    let targetSDF = engines[0].SDF;

    let data = [];
    engines.forEach(e => {
        const sdfDiff = e.compareSDF(targetSDF);
        const curvature = e.calculateCurvature();

        data.push({
            "SDF deviance": sdfDiff, 
            "Curve length": curvature[0], 
            "1st Deg Smoothness":  curvature[1] / curvature[0],
            "2nd Deg Smoothness":  curvature[2] / curvature[0]
        });
    });
    data.push([])

    console.table(data)
});
canvasA.addEventListener('pointermove', event => {
    const dimensions = canvasA.getBoundingClientRect();
    let mouseX = event.clientX - dimensions.left;
    let mouseY = event.clientY - dimensions.top;

    if(mouseDown){
        engines.map((e, i) => {
            e.updateStroke([mouseX, mouseY])
        });
    }
});

const cornerTest = () => {
    const correctPointSet = [];
    let w = engines[0].canvas.width;
    let h = engines[0].canvas.height;
    let points = [[0.3, 0.1], [0.3, 0.8], [0.8, 0.8]];
    for(let i = 0; i < points.length - 1; i++){
        let lerp = (a, b, t) => [
            w * (a[0] + t * (b[0] - a[0])) | 0,
            h * (a[1] + t * (b[1] - a[1])) | 0
        ];
    
        for(let t = 0; t <= 1; t += 0.01){
            correctPointSet.push(lerp(points[i], points[i + 1], t))
        }
    }
    engines.forEach(e => {
        e.createStroke();
        correctPointSet.forEach(point => e.updateStroke(point));
    });
    engines.forEach(e => {
        e.calculateSDF();
    });
    let targetSDF = engines[0].SDF;
    let data = [];
    engines.forEach(e => {
        const sdfDiff = e.compareSDF(targetSDF);
        const curvature = e.calculateCurvature();
    
        data.push({
            "SDF deviance": sdfDiff, 
            "Curve length": curvature[0], 
            "1st Deg Smoothness":  curvature[1] / curvature[0],
            "2nd Deg Smoothness":  curvature[2] / curvature[0]
        });
    });
    data.push([])
    console.table(data)
    console.log(data.length)
}
// cornerTest();

window.globalThis.cornerTest = cornerTest;

window.globalThis.saveRecentStroke = () => {
    return JSON.stringify(engines[3].finalPixels[engines[3].finalPixels.length - 1]);
}

const rawReferenceStroke = "[[85,157],[85,156],[85,156],[85,156],[85,156],[85,155],[85,155],[85,155],[85,154],[85,154],[85,153],[85,152],[85,152],[85,151],[85,151],[85,150],[85,150],[85,149],[85,148],[85,148],[85,147],[85,147],[86,146],[86,146],[86,145],[86,144],[86,144],[87,143],[87,143],[87,142],[88,141],[88,141],[89,141],[90,140],[90,140],[91,140],[92,139],[92,139],[93,138],[94,138],[95,137],[95,137],[96,137],[97,136],[97,136],[98,136],[99,135],[100,135],[100,135],[101,135],[102,134],[103,134],[103,134],[104,134],[105,134],[106,133],[107,133],[107,133],[108,133],[109,133],[110,133],[110,133],[111,133],[112,132],[112,132],[113,132],[114,132],[115,132],[115,132],[116,132],[117,132],[118,132],[118,132],[119,132],[120,132],[120,132],[121,132],[122,132],[123,132],[123,132],[124,132],[125,132],[125,132],[126,132],[127,132],[127,132],[128,132],[129,133],[130,133],[130,133],[131,133],[132,133],[133,133],[133,133],[134,133],[135,134],[135,134],[136,134],[137,135],[138,135],[138,135],[139,135],[140,136],[141,136],[141,136],[142,136],[143,137],[144,137],[144,137],[145,137],[146,138],[147,138],[147,138],[148,138],[149,138],[150,138],[150,138],[151,138],[152,139],[153,139],[153,139],[154,139],[155,140],[156,140],[156,140],[157,140],[158,140],[159,141],[160,141],[160,141],[161,141],[162,142],[163,142],[163,142],[164,142],[165,142],[166,143],[167,143],[167,143],[168,143],[169,143],[170,144],[171,144],[171,144],[172,144],[173,144],[174,145],[175,145],[175,145],[176,145],[177,145],[178,146],[179,146],[179,146],[180,146],[181,146],[182,147],[183,147],[184,147],[184,147],[185,147],[186,147],[187,147],[188,147],[189,147],[189,147],[190,147],[191,147],[192,148],[193,148],[193,148],[194,148],[195,148],[196,148],[197,148],[198,148],[198,148],[199,148],[200,148],[201,148],[202,148],[202,148],[203,148],[204,148],[205,149],[206,149],[206,149],[207,149],[208,149],[209,149],[210,149],[210,149],[211,149],[212,149],[213,149],[214,149],[214,149],[215,149],[216,149],[217,148],[218,148],[218,148],[219,148],[220,148],[221,148],[221,148],[222,148],[223,147],[224,147],[224,147],[225,147],[226,147],[227,147],[227,147],[228,147],[229,146],[230,146],[230,146],[231,146],[232,145],[232,145],[233,145],[234,144],[234,144],[235,144],[236,143],[236,143],[237,142],[238,141],[238,141],[239,141],[240,140],[240,140],[241,139],[241,139],[241,138],[242,137],[242,137],[242,136],[242,136],[243,135],[243,135],[243,134],[243,134],[242,133],[242,133],[242,132],[242,132],[241,131],[241,131],[240,130],[240,130],[239,130],[238,129],[238,129],[237,129],[236,129],[235,129],[235,129],[234,129],[233,128],[233,128],[232,128],[231,128],[230,128],[230,128],[229,128],[228,129],[228,129],[227,129],[226,129],[226,129],[225,129],[224,130],[223,130],[223,130],[222,130],[221,131],[221,131],[220,132],[219,133],[219,133],[218,133],[217,134],[217,134],[216,135],[215,136],[215,136],[214,137],[213,138],[212,139],[212,139],[211,140],[211,141],[210,142],[210,142],[210,143],[209,144],[209,145],[208,146],[208,146],[207,147],[207,148],[206,149],[205,150],[205,150],[204,151],[204,152],[203,153],[202,154],[202,154],[201,155],[201,156],[200,157],[200,158],[199,159],[199,159],[199,160],[198,161],[198,162],[197,163],[197,164],[196,165],[196,165],[195,166],[195,167],[194,168],[194,169],[193,170],[193,170],[192,171],[192,172],[191,173],[191,174],[190,175],[190,175],[190,176],[189,177],[189,178],[188,179],[188,180],[187,181],[187,181],[187,182],[186,183],[186,184],[185,185],[185,186],[184,187],[184,187],[183,188],[183,189],[182,190],[182,191],[181,192],[181,192],[180,193],[180,194],[179,195],[178,196],[178,197],[177,198],[177,198],[177,199],[176,200],[176,201],[175,202],[175,203],[174,204],[174,204],[173,205],[172,206],[172,207],[171,208],[170,209],[170,209],[170,210],[169,211],[169,212],[168,213],[168,214],[167,215],[167,215],[166,216],[165,217],[165,218],[164,219],[163,220],[163,220],[162,221],[161,222],[160,223],[159,224],[158,225],[158,225],[157,226],[156,227],[156,228],[155,229],[154,230],[154,230],[153,231],[152,232],[151,233],[150,234],[149,235],[149,235],[148,236],[147,237],[147,238],[146,239],[145,240],[145,240],[144,241],[143,242],[142,243],[141,244],[141,244],[140,245],[139,246],[138,246],[137,247],[136,248],[136,248],[135,249],[134,249],[133,250],[132,251],[132,251],[131,252],[130,253],[129,254],[129,254],[128,255],[127,255],[126,256],[125,257],[125,257],[124,257],[123,258],[122,258],[121,259],[121,259],[120,259],[119,260],[118,260],[117,261],[117,261],[116,261],[115,262],[114,262],[114,262],[113,262],[112,263],[111,263],[111,263],[110,263],[109,264],[108,264],[108,264],[107,264],[106,264],[105,264],[105,264],[104,264],[103,264],[103,264],[102,264],[101,264],[100,264],[100,264],[99,264],[99,264],[98,264],[97,264],[97,264],[96,264],[96,264],[95,263],[95,263],[94,263],[93,263],[93,263],[92,262],[92,262],[91,261],[91,261],[91,260],[91,259],[91,259],[90,258],[90,258],[90,257],[89,256],[89,256],[89,255],[89,254],[89,253],[89,253],[89,252],[89,251],[89,251],[89,250],[89,249],[89,249],[89,248],[90,247],[90,247],[90,246],[91,245],[91,244],[91,244],[92,243],[93,242],[93,242],[94,241],[95,240],[95,240],[96,239],[96,238],[97,237],[97,237],[98,236],[99,236],[100,235],[100,235],[101,235],[102,234],[103,234],[103,234],[104,234],[105,233],[106,233],[107,232],[107,232],[108,232],[109,231],[110,231],[110,231],[111,231],[112,231],[113,230],[114,230],[114,230],[115,230],[116,230],[117,230],[118,230],[118,230],[119,230],[120,229],[121,229],[121,229],[122,229],[123,229],[124,229],[125,229],[125,229],[126,229],[127,229],[128,230],[129,230],[130,230],[130,230],[131,230],[132,230],[133,230],[134,230],[134,230],[135,230],[136,230],[137,231],[138,231],[139,231],[139,231],[140,231],[141,232],[142,232],[143,233],[144,233],[144,233],[145,233],[146,233],[147,234],[148,234],[149,234],[149,234],[150,234],[151,235],[152,235],[153,236],[153,236],[154,236],[155,236],[156,237],[157,237],[158,237],[158,237],[159,237],[160,238],[161,238],[162,239],[162,239],[163,239],[164,240],[165,240],[166,241],[167,241],[167,241],[168,241],[169,242],[170,242],[171,243],[171,243],[172,243],[173,243],[174,244],[175,244],[175,244],[176,244],[177,245],[178,245],[179,246],[179,246],[180,246],[181,247],[182,247],[183,248],[183,248],[184,249],[185,249],[186,250],[186,250],[187,250],[188,250],[189,251],[190,251],[190,251],[191,251],[192,252],[193,252],[194,253],[194,253],[195,253],[196,254],[197,254],[197,254],[198,254],[199,255],[200,255],[201,256],[201,256],[202,256],[203,257],[204,257],[204,257],[205,257],[206,258],[207,258],[207,258],[208,258],[209,259],[210,259],[210,259],[211,259],[212,259],[213,260],[214,260],[214,260],[215,260],[216,261],[217,261],[217,261],[218,261],[219,262],[220,262],[220,262],[221,262],[222,262],[223,262],[224,262],[224,262],[225,262],[226,263],[227,263],[227,263],[228,263],[229,263],[230,263],[230,263],[231,263],[232,263],[233,263],[234,263],[234,263],[235,263],[236,263],[237,263],[238,263],[238,263],[239,263],[240,263],[241,263],[241,263],[242,263],[243,263],[244,263],[245,263],[245,263],[246,263],[247,263],[248,263],[248,263],[249,263],[250,262],[251,262],[251,262],[252,262],[253,261],[254,261],[254,261],[255,261],[256,261],[257,261],[257,261],[258,261],[259,260],[259,260],[260,260],[261,259],[261,259],[262,258],[263,257],[263,257],[264,256]]";
const referenceStroke = JSON.parse(rawReferenceStroke);
engines.map(e => e.createStroke())
engines.map(e => e.drawStrokes(referenceStroke))
engines.map(e => e.compareSDFV2(referenceStroke))