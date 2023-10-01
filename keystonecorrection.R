# Linear keystone distortion correction
# www.overfitting.net
# https://www.overfitting.net/2023/09/transformacion-trapezoidal-de-imagenes.html


# Keystone correction equations:
# https://discorpy.readthedocs.io/en/latest/tutorials/methods.html
# "2.2.7. Calculating coefficients of a correction model for perspective distortion"

# Distorted points (source)
xd=c(0, 900, 850, 200)
yd=c(0, -50, 650, 800)

# Undistorted points (destination)
# Example 1 (easy)
xu=c(200, 900, 950, 50) + 1100
yu=c(0, 200, 800, 650) + 500

# Example 2 (difficult)
xu=c(50, 100, 900, 800)+1100
yu=c(-400, 800, 600, -100)+200


# Solve 8 equations linear system: A * k = b -> k = inv(A) * b
A=matrix(nrow=8, ncol=8)
A[1,]=c(xd[1], yd[1], 1, 0,     0,     0, -xd[1]*xu[1], -yd[1]*xu[1])
A[2,]=c(0,     0,     0, xd[1], yd[1], 1, -xd[1]*yu[1], -yd[1]*yu[1])
A[3,]=c(xd[2], yd[2], 1, 0,     0,     0, -xd[2]*xu[2], -yd[2]*xu[2])
A[4,]=c(0,     0,     0, xd[2], yd[2], 1, -xd[2]*yu[2], -yd[2]*yu[2])
A[5,]=c(xd[3], yd[3], 1, 0,     0,     0, -xd[3]*xu[3], -yd[3]*xu[3])
A[6,]=c(0,     0,     0, xd[3], yd[3], 1, -xd[3]*yu[3], -yd[3]*yu[3])
A[7,]=c(xd[4], yd[4], 1, 0,     0,     0, -xd[4]*xu[4], -yd[4]*xu[4])
A[8,]=c(0,     0,     0, xd[4], yd[4], 1, -xd[4]*yu[4], -yd[4]*yu[4])

b=as.matrix(c(xu[1], yu[1], xu[2], yu[2], xu[3], yu[3], xu[4], yu[4]))

k=solve(A, b)  # equivalent to inv(A) * b = solve(A) %*% b

# Undo distortion function
undo.keystone = function(xd, yd, k) {
    xu = (k[1]*xd+k[2]*yd+k[3]) / (k[7]*xd+k[8]*yd+1)
    yu = (k[4]*xd+k[5]*yd+k[6]) / (k[7]*xd+k[8]*yd+1)
    return(c(xu, yu))  # return pair (xu, yu)
}

# Check
for (i in 1:4) {
    print(undo.keystone(xd[i], yd[i], k))
}


# Plot trapezoids
plot(c(xd, xd[1]), c(yd, yd[1]), type='l', col='blue', asp=1,
     xlab='X', ylab='Y', xlim=c(0,2100), ylim=c(0,800))
lines(c(xu, xu[1]), c(yu, yu[1]), type='l', col='red')
for (i in 1:4) {
    lines(c(xd[i], xu[i]), c(yd[i], yu[i]), type='l', lty=3, col='darkgray')
}


# Test point(s)
xdp=c(180)
ydp=c(600)

lines(xdp, ydp, type='p', col='blue', pch=16, cex=1.5)
for (i in 1:length(xdp)) {
    xuyu=undo.keystone(xdp[i], ydp[i], k)
    lines(xuyu[1], xuyu[2], type='p', col='red', pch=16, cex=1.5)
}
