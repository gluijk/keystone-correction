# Keystone distortion correction
# www.overfitting.net
# https://www.overfitting.net/2025/04/transformacion-trapezoidal-de-imagenes.html

library(tiff)


# Keystone correction equations:
# https://discorpy.readthedocs.io/en/latest/tutorials/methods.html
# "2.2.7. Calculating coefficients of a correction model for perspective distortion"

# Distorted points (source)
imgd=readTIFF("distorted.tif")
xu=c(266, 1132, 3470, 1869)  # top-left, bottom-left, bottom-right, top-right
yu=c(1130, 3186, 1719, 692)

# Undistorted points (destination)
posx=(xu[1]+xu[2])/2
posy=(yu[1]+yu[4])/2

# Real aspect ratio of Football field:
halflength=378.3
width=454.6

f=4.7  # scaling factor
xd=c(posx, posx,         posx+halflength*f, posx+halflength*f)
yd=c(posy, posy+width*f, posy+width*f,      posy)

# NOTE: we swap the distorted and undistorted trapezoids because
# we want to model the transformation
# FROM CORRECTED coords (DST) -> TO UNCORRECTED coords (ORG)


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
    xu=(k[1]*xd+k[2]*yd+k[3]) / (k[7]*xd+k[8]*yd+1)
    yu=(k[4]*xd+k[5]*yd+k[6]) / (k[7]*xd+k[8]*yd+1)
    return(c(xu, yu))  # return pair (xu, yu)
}

# Check
for (i in 1:4) {
    print(undo.keystone(xd[i], yd[i], k))
}


# Plot trapezoids
plot(c(xd, xd[1]), c(yd, yd[1]), type='l', col='red', asp=1,
     xlab='X', ylab='Y', xlim=c(1, ncol(imgd)), ylim=c(nrow(imgd), 1))
lines(c(xu, xu[1]), c(yu, yu[1]), type='l', col='blue')
for (i in 1:4) {
    lines(c(xd[i], xu[i]), c(yd[i], yu[i]), type='l', lty=3, col='darkgray')
}
abline(h=c(1,nrow(imgd)), v=c(1,ncol(imgd)))


# Correct keystone distortion
DIMXd=ncol(imgd)
DIMYd=nrow(imgd)

imgc=imgd*0
DIMXc=ncol(imgc)
DIMYc=nrow(imgc)
for (x in 1:DIMXc) {
    for (y in 1:DIMYc) {
        xuyu=round(undo.keystone(x, y, k))
        if (xuyu[1]>=1 & xuyu[1]<=DIMXd & xuyu[2]>=1 & xuyu[2]<=DIMYd)
            imgc[y, x,]=imgd[xuyu[2], xuyu[1],]  # nearest neighbour interp
    }
}
writeTIFF(imgc, "corrected.tif", bits.per.sample=16)
