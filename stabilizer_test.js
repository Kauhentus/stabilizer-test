const canvasA = document.getElementById('canvasA');
const canvasB = document.getElementById('canvasB');
const canvasC = document.getElementById('canvasC');
const canvasD = document.getElementById('canvasD');
const canvasE = document.getElementById('canvasE');
const canvasF = document.getElementById('canvasF');

const stabilizerConstructor = (processPointFunction, id) => {
    return {
        points: [],
        addPoint: (point) => {
            stabilizers[id].points.push(point);
        },
        returnPoints: () => {
            return processPointFunction(stabilizers[id].points);
        },
        reset: () => {
            stabilizers[id].points = [];
        }
    };
}

const stabilizers = [
    // identity
    stabilizerConstructor((points) => {
        return points;
    }, 0),

    // interpolation
    stabilizerConstructor((points) => {
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
    }, 1),

    /*stabilizerConstructor((points) => {
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

        return collectedPoints;
    }),*/

    // bezier 
    stabilizerConstructor((points) => {
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
    }, 2),

    // calculate vel & acc
    /*stabilizerConstructor((points) => {
        let prevPoint = points[0];
        let velocities = [];

        for(let i = 0; i < points.length; i++){
            let currentPoint = points[i];

            let velocity = [
                (currentPoint[0] - prevPoint[0]) / 1,
                (currentPoint[1] - prevPoint[1]) / 1,
            ];
            velocities.push(velocity);

            prevPoint = currentPoint;
        }

        let prevVelocity = velocities[0];
        let accelerations = [];
        for(let i = 0; i < points.length; i++){
            let currentVelocity = velocities[i];

            let acceleration = [
                (currentVelocity[0] - prevVelocity[0]) / 1,
                (currentVelocity[1] - prevVelocity[1]) / 1,
            ];
            accelerations.push(acceleration);

            prevVelocity = currentVelocity;
        }

        let finalPoints = points;

        return finalPoints;
    }, 3),*/

    // spring model
    stabilizerConstructor((points) => {
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
    }, 3),

    // savitsky-golay 5th-degree
    stabilizerConstructor((points) => {

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
    }, 4),

    // savitsky-golay 9th-degree
    stabilizerConstructor((points) => {

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
    }, 5),
    
    // catmull rom
    /*stabilizerConstructor((points) => {
        const maxDistance = 50;
        const calcPointDist = (p1, p2) => Math.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2);
    
        let collectedPoints = [
            [
                points[0][0] + 1,
                points[0][1] + 1,
            ], 
            points[0]
        ];
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
        collectedPoints.push([actualLast[0] + 1, actualLast[1] + 1]);

        const alpha = 0.5;
        const tension = 0; 

        const calcSegment = (p0, p1, p2, p3) => {
            const getDist = (a, b) => Math.sqrt((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2);
            const addVec = (a, b) => [a[0] + b[0], a[1] + b[1]];
            const subVec = (a, b) => [a[0] - b[0], a[1] - b[1]];
            const mulVec = (a, s) => [a[0] * s, a[1] * s];
            const divVec = (a, s) => [a[0] / s, a[1] / s];

            const t0 = 0;
            const t1 = t0 + getDist(p0, p1) ** alpha; 
            const t2 = t1 + getDist(p1, p2) ** alpha; 
            const t3 = t2 + getDist(p2, p3) ** alpha; 

            const m1 = mulVec(
                addVec(
                    subVec(
                        divVec(subVec(p1, p0), t1 - t0), 
                        divVec(subVec(p2, p0), t2 - t0)
                    ),
                    divVec(subVec(p2, p1), t2 - t1)
                ),
                (1 - tension) * (t2 - t1)
            );

            const m2 = mulVec(
                addVec(
                    subVec(
                        divVec(subVec(p2, p1), t2 - t1), 
                        divVec(subVec(p3, p1), t3 - t1)
                    ),
                    divVec(subVec(p3, p2), t3 - t2)
                ),
                (1 - tension) * (t2 - t1)
            );

            const a = addVec(
                addVec(
                    mulVec(subVec(p1, p2), 2),
                    m1
                ),
                m2
            );

            const b = subVec(
                subVec(
                    subVec(
                        mulVec(subVec(p1, p2), -3),
                        m1
                    ),
                    m1
                ),
                m2
            );

            const c = m1;
            const d = p1;

            let points = [];
            for(let t = 0; t <= 1; t += 0.2){
                points.push(addVec(
                    addVec(
                        addVec(
                            mulVec(a, t * t * t),
                            mulVec(b, t * t)
                        ),
                        mulVec(c, t)
                    ),
                    d
                ));
            }

            return points;
        }

        const calcChain = (points = []) => {
            if(points.length < 4) return points;

            points.unshift(points[0]);
            points.push(points[points.length - 1]);

            let newPoints = [];
            for(let i = 0; i < points.length - 3; i++){
                let curPoints = calcSegment(
                    points[i], 
                    points[i + 1], 
                    points[i + 2], 
                    points[i + 3], 
                    0.5
                );

                newPoints.push(...curPoints);
            };

            return newPoints;
        };

        return calcChain(collectedPoints);
    })*/
];

let mouseDown = false;
canvasA.addEventListener('pointerdown', () => mouseDown = true);
canvasA.addEventListener('pointerup', () => {
    mouseDown = false;
    buffer = [];
    stabilizers.map(s => s.reset());
});

let buffer = [];
let mouseX = 0, mouseY = 0;
canvasA.addEventListener('pointermove', event => {
    const dimensions = canvasA.getBoundingClientRect();
    mouseX = event.clientX - dimensions.left;
    mouseY = event.clientY - dimensions.top;

    console.log("REE")

    if(mouseDown){
        // console.log(mouseX, mouseY)
        buffer.push([mouseX, mouseY]);
        broadcast([mouseX, mouseY]);
    }
});

const canvasContainers = [canvasA, canvasB, canvasC, canvasD, canvasE, canvasF];

const canvasObjects = canvasContainers.map(canvasContainer => {
    const canvasDimensions = canvasContainer.getBoundingClientRect();
    const width = canvasDimensions.width;
    const height = canvasDimensions.height;

    const canvasObject = document.createElement('canvas');
    canvasObject.width = width;
    canvasObject.height = height;
    canvasContainer.appendChild(canvasObject);

    const canvasContext = canvasObject.getContext('2d');

    return {
        canvas: canvasObject,
        ctx: canvasContext
    };
});

const engines = canvasObjects.map(pair => {
    let canvas = pair.canvas;
    let ctx = pair. ctx;
    let lastPoint = -1;

    return {
        addPoint: point => {
            // ctx.fillRect(point[0], point[1], 2, 2);
            if(lastPoint == -1) lastPoint = point;

            ctx.beginPath();
            ctx.moveTo(lastPoint[0], lastPoint[1]);
            ctx.lineTo(point[0], point[1]);
            ctx.closePath();
            ctx.stroke();

            lastPoint = point;
        },
        clear: () => {
            ctx.clearRect(0, 0, canvas.width, canvas.height);
        }
    }
});

const broadcast = (point) => {
    engines.forEach((engine, i) => {
        const currentStabilizer = stabilizers[i];
        currentStabilizer.addPoint(point);
        const currentPoints = currentStabilizer.returnPoints();

        engine.clear();

        currentPoints.forEach(p => engine.addPoint(p));
    });
};

// gaussian elimination from https://rosettacode.org/wiki/Gaussian_elimination#JavaScript
// Lower Upper Solver
function lusolve(A, b, update) {
	var lu = ludcmp(A, update)
	if (lu === undefined) return // Singular Matrix!
	return lubksb(lu, b, update)
}
 
// Lower Upper Decomposition
function ludcmp(A, update) {
	// A is a matrix that we want to decompose into Lower and Upper matrices.
	var d = true
	var n = A.length
	var idx = new Array(n) // Output vector with row permutations from partial pivoting
	var vv = new Array(n)  // Scaling information
 
	for (var i=0; i<n; i++) {
		var max = 0
		for (var j=0; j<n; j++) {
			var temp = Math.abs(A[i][j])
			if (temp > max) max = temp
		}
		if (max == 0) return // Singular Matrix!
		vv[i] = 1 / max // Scaling
	}
 
	if (!update) { // make a copy of A 
		var Acpy = new Array(n)
		for (var i=0; i<n; i++) {		
			var Ai = A[i] 
			Acpyi = new Array(Ai.length)
			for (j=0; j<Ai.length; j+=1) Acpyi[j] = Ai[j]
			Acpy[i] = Acpyi
		}
		A = Acpy
	}
 
	var tiny = 1e-20 // in case pivot element is zero
	for (var i=0; ; i++) {
		for (var j=0; j<i; j++) {
			var sum = A[j][i]
			for (var k=0; k<j; k++) sum -= A[j][k] * A[k][i];
			A[j][i] = sum
		}
		var jmax = 0
		var max = 0;
		for (var j=i; j<n; j++) {
			var sum = A[j][i]
			for (var k=0; k<i; k++) sum -= A[j][k] * A[k][i];
			A[j][i] = sum
			var temp = vv[j] * Math.abs(sum)
			if (temp >= max) {
				max = temp
				jmax = j
			}
		}
		if (i <= jmax) {
			for (var j=0; j<n; j++) {
				var temp = A[jmax][j]
				A[jmax][j] = A[i][j]
				A[i][j] = temp
			}
			d = !d;
			vv[jmax] = vv[i]
		}
		idx[i] = jmax;
		if (i == n-1) break;
		var temp = A[i][i]
		if (temp == 0) A[i][i] = temp = tiny
		temp = 1 / temp
		for (var j=i+1; j<n; j++) A[j][i] *= temp
	}
	return {A:A, idx:idx, d:d}
}
 
// Lower Upper Back Substitution
function lubksb(lu, b, update) {
	// solves the set of n linear equations A*x = b.
	// lu is the object containing A, idx and d as determined by the routine ludcmp.
	var A = lu.A
	var idx = lu.idx
	var n = idx.length
 
	if (!update) { // make a copy of b
		var bcpy = new Array(n) 
		for (var i=0; i<b.length; i+=1) bcpy[i] = b[i]
		b = bcpy
	}
 
	for (var ii=-1, i=0; i<n; i++) {
		var ix = idx[i]
		var sum = b[ix]
		b[ix] = b[i]
		if (ii > -1)
			for (var j=ii; j<i; j++) sum -= A[i][j] * b[j]
		else if (sum)
			ii = i
		b[i] = sum
	}
	for (var i=n-1; i>=0; i--) {
		var sum = b[i]
		for (var j=i+1; j<n; j++) sum -= A[i][j] * b[j]
		b[i] = sum / A[i][i]
	}
	return b // solution vector x
}