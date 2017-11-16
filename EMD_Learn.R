
###############
# EMD（Empirical Mode Decomposition, 经验模态分解） 学习
# R Package : EMD
# JIANGWEI 
# 2017-11-15



library(EMD)


# 1、cvtype ------------------------------------------------------------------

# cvtype Generating test dataset index for cross-validation 
# 返回的矩阵每行是一个测试集的索引

# Traditional 4-fold cross-validation for 100 observations
cvtype(n=100, cv.bsize=1, cv.kfold=4, cv.random=FALSE)
# Random 4-fold cross-validation with block size 2 for 100 observations
cvtype(n=100, cv.bsize=2, cv.kfold=4, cv.random=TRUE)




# 2、emd ---------------------------------------------------------------------

# emd  经验模态分解

### Empirical Mode Decomposition
ndata <- 3000
tt2 <- seq(0, 9, length=ndata)
xt2 <- sin(pi * tt2) + sin(2* pi * tt2) + sin(6 * pi * tt2)  + 0.5 * tt2

try <- emd(xt2, tt2, boundary="wave",plot.imf=FALSE)

### Ploting the IMF's
par(mfrow=c(try$nimf+2, 1), mar=c(2,1,2,1))
rangeimf <- range(try$imf)
for(i in 1:try$nimf) {
  plot(tt2, try$imf[,i], type="l", xlab="", ylab="", ylim=rangeimf,
       main=paste(i, "-th IMF", sep="")); abline(h=0)
}
plot(tt2, try$residue, xlab="", ylab="", main="residue", type="l", axes=TRUE);box()
plot(tt2, sin(2*pi * tt2), xlab="", ylab="", main="对比：sin(2*pi * tt2)", type="l", axes=T);box() 




# 3、emd.pred --------------------------------------------------------------

# This function calculates prediction values and confidence limits using EMD and VAR (vector autoregressive)model.
# Usage
# emd.pred(varpred, trendpred, ci = 0.95, figure = TRUE)
# 
# Arguments
# varpred    prediction result of IMF’s by VAR model.
# trendpred  prediction result of residue by polynomial regression model.
# ci         confidence interval level.
# figure     specifies whether prediction result is displayed.




# 4、emd2d -----------------------------------------------------------------

# This function performs the bidimenasional empirical mode decomposition utilizing extrema detection
# based on the equivalence relation between neighboring pixels.
# 该函数以相邻像素之间的等价关系为基础，利用极值检测方法进行二维经验模式分解。

data(lena)
z <- lena[ seq(1, 512, by=16), seq(1, 512, by=16) ]
image(lena, main="原始图像lena", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)
image(z, main="约减后的图像lena", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)

#下面是二维经验模态分解函数调用，可能需运行较长时间
lenadecom <- emd2d(z, max.imf = 4)
imageEMD(z=z, emdz=lenadecom, extrema=TRUE, col=gray(0:100/100))


### 例二：Test Image
ndata <- 128

#模拟产生两个和余额减后图像大小相同的数据集meanf1、meanf2以及它们的和meanf
x <- y <- seq(0, 9, length=ndata)
meanf1 <- outer(sin(2 * pi * x), sin(2 * pi * y))
meanf2 <- outer(sin(0.5 * pi * x), sin(0.5 * pi * y))
meanf <- meanf1 + meanf2

#又模拟多产生一个噪音数据，加上去
snr <- 2
set.seed(77)
zn <- meanf + matrix(rnorm(ndata^2, 0, sd(c(meanf))/snr), ncol=ndata)


rangezn <- range(c(meanf1, meanf2, meanf, zn)) #上述数据的区间

#分别把上面模拟产生的四个数据集绘制成图像
par(mfrow=c(2,2), mar=0.1 + c(0, 0.25, 3, 0.25))
image(meanf1, main="high frequency component", xlab="", ylab="", zlim=rangezn, 
      col=gray(100:0/100), axes=FALSE)
image(meanf2, main="low frequency component", xlab="", ylab="", zlim=rangezn, 
      col=gray(100:0/100), axes=FALSE)
image(meanf, main="test image", xlab="", ylab="", zlim=rangezn, col=gray(100:0/100), axes=FALSE)
image(zn, main="noisy image", xlab="", ylab="", zlim=rangezn, col=gray(100:0/100), axes=FALSE)

#下面是二维经验模态分解函数调用，可能需运行较长时间
out <- emd2d(zn, max.imf=3, sm="locfit", smlevels=1, spar=0.004125)
par(mfcol=c(3,1), mar=0.1 + c(0, 0.25, 0.25, 0.25))
image(out$imf[[1]], main="", xlab="", ylab="", col=gray(100:0/100), zlim=rangezn, axes=FALSE)
image(out$imf[[2]], main="", xlab="", ylab="", col=gray(100:0/100), zlim=rangezn, axes=FALSE)
image(out$imf[[3]], main="", xlab="", ylab="", col=gray(100:0/100), zlim=rangezn, axes=FALSE)




# 5、emddenoise ------------------------------------------------------------

# 该函数通过经验模式分解和阈值进行去噪。

ndata <- 1024
tt <- seq(0, 9, length=ndata)
#三段不同频率的数据叠加
meanf <- (sin(pi*tt) + sin(2*pi*tt) + sin(6*pi*tt)) * (0.0<tt & tt<=3.0) +
  (sin(pi*tt) + sin(6*pi*tt)) * (3.0<tt & tt<=6.0) +
  (sin(pi*tt) + sin(6*pi*tt) + sin(12*pi*tt)) * (6.0<tt & tt<=9.0)

snr <- 3.0
#三个新的标准差
sigma <- c(sd(meanf[tt<=3]) / snr, sd(meanf[tt<=6 & tt>3]) / snr,
           sd(meanf[tt>6]) / snr)
set.seed(1)
#模拟新的三个噪音数据集
error <- c(rnorm(sum(tt<=3), 0, sigma[1]),
           rnorm(sum(tt<=6 & tt>3), 0, sigma[2]), 
           rnorm(sum(tt>6), 0, sigma[3]))
#合成数据：数据+噪音
xt <- meanf + error
#交叉验证索引，测试集要用到
cv.index <- cvtype(n=ndata, cv.kfold=2, cv.random=FALSE)$cv.index

try10 <- emddenoise(xt, cv.index=cv.index, cv.level=2, by.imf=TRUE)
try10$optlambda

# 返回的结果包含以下信息
# dxt          : denoised signal
# optlambda    : threshold values by cross-validation
# lambdaconv   : sequence of lambda’s by cross-validation
# perr         : sequence of prediction error by cross-validation
# demd         : denoised IMF’s and residue
# niter        : the number of iteration for optimal threshold value





# 6、extractimf ------------------------------------------------------------

# 该函数从给定信号提取固有模态函数（基本模式分量imf）。
# 给定信号可以是原始信号、残差等
# （在希尔伯特黄变换中本征模态函数是基于序列数据的局部时间尺度特征而得出）

### Generating a signal
ndata <- 3000
par(mfrow=c(1,1), mar=c(1,1,1,1))
tt2 <- seq(0, 9, length=ndata)
xt2 <- sin(pi * tt2) + sin(2* pi * tt2) + sin(6 * pi * tt2) + 0.5 * tt2
plot(tt2, xt2, xlab="", ylab="", type="l", axes=FALSE); box()
### Extracting the first IMF by sifting process
tryimf <- extractimf(xt2, tt2, check=FALSE)  #提取第一分量




# 7、extractimf2d ----------------------------------------------------------

# This function extracts the bidimensional intrinsic mode function from given an image utilizing
# extrema detection based on the equivalence relation between neighboring pixels.
# 该函数根据相邻像素之间的等价关系，利用极值检测方法提取二维固有模态函数。

data(lena)
z <- lena[seq(1, 512, by=16), seq(1, 512, by=16)]
# 提取第一分量，是emd2d的更方便的一个变种函数。emd2d是给出来所有的分量
lenaimf1 <- extractimf2d(z, check=FALSE) 




# 8、extrema ---------------------------------------------------------------
# This function indentifies extrema and zero-crossings.
# 这个函数可以识别局部极值点和零交叉点的位置。
# 返回结果的形式是矩阵：每个位置用两个点的索引表示，如，交叉的话，则用前后两个点的索引表示，
# 若刚好在零线上，则两个该点的索引

y <- c(0, 1, 2, 1, -1, 1:4, 5, 6, 0, -4, -6, -5:5, -2:2)
#y <- c(0, 0, 0, 1, -1, 1:4, 4, 4, 0, 0, 0, -5:5, -2:2, 2, 2)
#y <- c(0, 0, 0, 1, -1, 1:4, 4, 4, 0, 0, 0, -5:5, -2:2, 0, 0)
plot(y, type = "b"); abline(h = 0)
extrema(y)




# 9、extrema2dC ------------------------------------------------------------
# This function finds the bidimensional local extrema based on the equivalence relation between
# neighboring pixels.

data(lena)
z <- lena[seq(1, 512, by=4), seq(1, 512, by=4)]
par(mfrow=c(1,3), mar=c(0, 0.5, 2, 0.5))
#原始图像
image(z, main="Lena", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)

#找到矩阵极值点
example <- extrema2dC(z=z)
localmin <- matrix(256, nrow = 128, ncol = 128)
#把原来图像的值放到新的空白矩阵中
localmin[example$minindex] <- z[example$minindex] 
image(localmin, main="Local minimum", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)

localmax <- matrix(0, 128, 128)
#把原图像中的最大值放入底色为黑色的空白矩阵中
localmax[example$maxindex] <- z[example$maxindex]
image(localmax, main="Local maximum", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)




# 10、hilbertspec ----------------------------------------------------------
# 该函数利用Hilbert变换计算振幅和瞬时频率。
# 函数返回值：
# amplitude   ： matrix of amplitudes for multiple signals xt
# instantfreq ： matrix of instantaneous frequencies for multiple signals xt
# energy      ： cumulative energy of multiple signals

tt <- seq(0, 0.1, length = 2001)[1:2000]
f1 <- 1776; f2 <- 1000
# 两个数据合成，第一个分布在两头，第二个在中间
xt <- sin(2*pi*f1*tt) * (tt <= 0.033 | tt >= 0.067) + sin(2*pi*f2*tt)

### Before treating intermittence
interm1 <- emd(xt, tt, boundary="wave", max.imf=2, plot.imf=FALSE)
### After treating intermittence
# interm=0.0007：指定将被排除在IMF之外的周期向量，以应对模式混合。
interm2 <- emd(xt, tt, boundary="wave", max.imf=2, plot.imf=FALSE, interm=0.0007)

par(mfrow=c(2,1), mar=c(2,2,2,1))
test1 <- hilbertspec(interm1$imf)
spectrogram(test1$amplitude[,1], test1$instantfreq[,1]) #x轴振幅，y轴频率
test2 <- hilbertspec(interm2$imf, tt=tt)
spectrogram(test2$amplitude[,1], test2$instantfreq[,1])



# 11、imageEMD -------------------------------------------------------------

# This function draws plots of input image, IMF’s, residue and extrema.
# 这个函数把原始图像、分量、残差、极值都画出来，方便对比

data(lena)
z <- lena[seq(1, 512, by=16), seq(1, 512, by=16)]
image(z, main="Lena", xlab="", ylab="", col=gray(0:100/100), axes=FALSE)
lenadecom <- emd2d(z, max.imf = 4)
imageEMD(z=z, emdz=lenadecom, extrema=TRUE, col=gray(0:100/100))



# 12、semd 统计经验模态分解 --------------------------------------------------------

# This function performs empirical mode decomposition using spline smoothing not interpolation for sifting process. 
# The smoothing parameter is automatically detemined by cross-validation.
# 该函数采用样条平滑法进行了经验模态分解，而不是对筛选过程进行插值。平滑参数被交叉验证自动检测。

ndata <- 2048
tt <- seq(0, 9, length=ndata)
# 合成信号
xt <- sin(pi * tt) + sin(2* pi * tt) + sin(6 * pi * tt) + 0.5 * tt
set.seed(1)
# 添加噪音
xt <- xt + rnorm(ndata, 0, sd(xt)/5)

#两种不同的光滑方法
### Empirical Mode Decomposition by Interpolation
emdbyint <- emd(xt, tt, max.imf = 5, boundary = "wave")
### Empirical Mode Decomposition by Smoothing
emdbysm <- semd(xt, tt, cv.kfold=4, boundary="wave", smlevels=1, max.imf=5)

#第一种方法（插值）结果绘图
par(mfcol=c(6,2), mar=c(2,2,2,1), oma=c(0,0,2,0))
rangext <- range(xt); rangeimf <- rangext - mean(rangext)
plot(tt, xt, xlab="", ylab="", main="signal", ylim=rangext, type="l")
mtext("Decomposition by EMD", side = 3, line = 2, cex=0.85, font=2)
plot(tt, emdbyint$imf[,1], xlab="", ylab="", main="imf 1", ylim=rangeimf, type="l")
abline(h=0, lty=2)
plot(tt, emdbyint$imf[,2], xlab="", ylab="", main="imf 2", ylim=rangeimf, type="l")
abline(h=0, lty=2)
plot(tt, emdbyint$imf[,3], xlab="", ylab="", main="imf 3", ylim=rangeimf, type="l")
abline(h=0, lty=2)
plot(tt, emdbyint$imf[,4], xlab="", ylab="", main="imf 4", ylim=rangeimf, type="l")
abline(h=0, lty=2)
plot(tt, emdbyint$imf[,5]+emdbyint$residue, xlab="", ylab="", main="remaining signal",
ylim=rangext, type="l")

#第二种方法结果（样条平滑）绘图
plot(tt, xt, xlab="", ylab="", main="signal", ylim=rangext, type="l")
mtext("Decomposition by SEMD", side = 3, line = 2, cex=0.85, font=2)
plot(tt, emdbysm$imf[,1], xlab="", ylab="", main="noise", ylim=rangeimf, type="l")
abline(h=0, lty=2)
plot(tt, emdbysm$imf[,2], xlab="", ylab="", main="imf 1", ylim=rangeimf, type="l")
abline(h=0, lty=2)
plot(tt, emdbysm$imf[,3], xlab="", ylab="", main="imf 2", ylim=rangeimf, type="l")
abline(h=0, lty=2)
plot(tt, emdbysm$imf[,4], xlab="", ylab="", main="imf 3", ylim=rangeimf, type="l")
abline(h=0, lty=2)
plot(tt, emdbysm$residue, xlab="", ylab="", main="residue", ylim=rangext, type="l")



# 13、spectrogram ----------------------------------------------------------

# This function produces image of amplitude by time index and instantaneous frequency. The horizontal
# axis represents time, the vertical axis is instantaneous frequency, and the color of each point
# in the image represents amplitude of a particular frequency at a particular time.

tt <- seq(0, 0.1, length = 2001)[1:2000]
f1 <- 1776; f2 <- 1000
xt <- sin(2*pi*f1*tt) * (tt <= 0.033 | tt >= 0.067) + sin(2*pi*f2*tt)
### Before treating intermittence
interm1 <- emd(xt, tt, boundary="wave", max.imf=2, plot.imf=FALSE)
### After treating intermittence
interm2 <- emd(xt, tt, boundary="wave", max.imf=2, plot.imf=FALSE,
               interm=0.0007)
par(mfrow=c(2,1), mar=c(2,2,2,1))
test1 <- hilbertspec(interm1$imf)
spectrogram(test1$amplitude[,1], test1$instantfreq[,1])
test2 <- hilbertspec(interm2$imf, tt=tt)
spectrogram(test2$amplitude[,1], test2$instantfreq[,1])



# 14、kospi200 数据集 ---------------------------------------------------------

# the weekly KOSPI 200 index from January, 1990 to February, 2007.

data(kospi200)
names(kospi200)
plot(kospi200$date, kospi200$index, type="l")



# 15、lena 数据集 -------------------------------------------------------------

# Gray Lena image：  A 512x512 gray image of Lena.

data(lena)
image(lena, col=gray(0:100/100), axes=FALSE)



# 16、lennon 数据集 -----------------------------------------------------------

# Gray John Lennon image ：A 256x256 gray image of John Lennon.

data(lennon)
image(lennon, col=gray(100:0/100), axes=FALSE)



# 17、solar 太阳辐射照度数据集 ----------------------------------------------------

# 三个不同作者提供的不同的太阳照度数据集
# solar irradiance proxy data.
# Hoyt and Schatten (1993) reconstructed solar irradiance (from 1700 through 1997) using the amplitude
# of the 11-year solar cycle together with a long term trend estimated from solar-like stars. They
# put relatively more weight on the length of the 11-year cycle.
# Lean et al. (1995) reconstructed solar irradiance (from 1610 through 2000) using the amplitude of
# the 11-year solar cycle and a long term trend estimated from solar-like stars.
# 10-Beryllium (10Be) is measured in polar ice from 1424 through 1985. 10-Beryllium (10Be) is
# produced in the atmosphere by incoming cosmic ray flux, which in turn is influenced by the solar
# activity. The higher the solar activity, the lower the flux of cosmic radiation entering the earth
# atmosphere and therefore the lower the production rate of 10Be. The short atmospheric lifetime of
# 10Be of one to two years (Beer et al. 1994) allows the tracking of solar activity changes and offers
# an alternative way to the sunspot based techniques for the analysis of the amplitude and length of
# the solar cycle as well as for low frequency variations.


data(solar.hs)
names(solar.hs)
plot(solar.hs$year, solar.hs$solar, type="l")

data(solar.lean)
names(solar.lean)
plot(solar.lean$year, solar.lean$solar, type="l")

data(beryllium)
names(beryllium)
plot(beryllium$year, beryllium$be, type="l")



# 18、sunspot 数据集----------------------------------------------------------

# sunspot from 1610 through 1995.从1610年到1995年的太阳黑子。

data(sunspot)
names(sunspot)
plot(sunspot$year, sunspot$sunspot, type="l")
