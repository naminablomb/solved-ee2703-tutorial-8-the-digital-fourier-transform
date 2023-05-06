Download Link: https://assignmentchef.com/product/solved-ee2703-tutorial-8-the-digital-fourier-transform
<br>
You are all studying DSP this semester, and you know the following from this course and the signals course:

<ul>

 <li>A time function of the type <em>f</em>(<em>t</em>)<em>u</em><sub>0</sub>(<em>t</em>) that grows slower than <em>e<sup>αt </sup></em>for some <em>α </em>has a Laplace transform, <em>F</em>(<em>s</em>). The Laplace transform has an inverse transform which involves integrating vertically in the complex plane along a path that is to the right of all poles.</li>

 <li>A finite energy function of the type <em>f</em>(<em>t</em>) has a Fourier Transform (and its inverse):</li>

</ul>

<em>F</em><em><sup>t</sup>dt</em>

<em>f</em>

<ul>

 <li>If <em>f</em>(<em>t</em>) is periodic with period 2<em>π</em>, the Fourier Transform collapses to the Fourier Series</li>

</ul>

<em>f</em><em> c<sub>n</sub>e</em><em>jnt</em>

<em>n</em>=−∞

where

<em>c</em><em>n </em><em> f</em>(<em>t</em>)<em>e</em>−<em>jntdt</em>

for any <em>t</em><sub>0</sub>. We can take <em>t</em><sub>0 </sub>= 0 for convenience.

<ul>

 <li>We can invert this picture and say, suppose <em>f</em>[<em>n</em>] are the samples of some function <em>f</em>(<em>t</em>), then we define the <em>Z </em>transform as</li>

</ul>

<em>F</em><em>n</em>

<em>n</em>=−∞

Replacing <em>z </em>with <em>e<sup>jθ </sup></em>we get

∞

<em>F</em>(<em>e</em><em>jθ</em>) = ∑ <em>f</em>[<em>n</em>]<em>e</em>−<em>jnθ n</em>=−∞

So clearly <em>F</em>(<em>z</em>) is like the periodic time function that gives rise to the fourier series whose coefficients are the samples <em>f</em>[<em>n</em>].

<em>F</em>(<em>e<sup>jθ</sup></em>) is called the Digital Spectrum of the samples <em>f</em>[<em>n</em>]. It is also called the DTFT of <em>f</em>[<em>n</em>]

<ul>

 <li><em>F</em>(<em>e<sup>jθ</sup></em>) is continuous and periodic. <em>f</em>[<em>n</em>] is discrete and aperiodic. Suppose now <em>f</em>[<em>n</em>] is itself periodic with a period <em>N</em>, i.e.,</li>

</ul>

<em>f</em>[<em>n</em>+<em>N</em>] = <em>f</em>[<em>n</em>] ∀<em>n</em>

Then, it should have samples for its DTFT. This is true, and leads to the Discrete Fourier Transform or the DFT:

Suppose <em>f</em>[<em>n</em>] is a periodic sequence of samples, with a period <em>N</em>. Then the DTFT of the sequence is also a periodic sequence <em>F</em>[<em>k</em>] with the same period <em>N</em>. (You will prove this in the DSP course). So we have

−1                                         <em>nk            </em><em>N</em>−1                  <em>nk</em>

<em>Ff</em>[<em>n</em>]<em>W</em>

<em>n</em>=0                                           <em>N            </em><em>n</em>=0

<em>N</em>−1

1      −<em>nk f</em>[<em>n</em>]   = ∑ <em>F</em>[<em>k</em>]<em>W</em>

<em>N </em><em>k</em>=0

Here <em>W </em> is used simply to make the equations less cluttered.

The values <em>F</em>[<em>k</em>] are what remains of the Digital Spectrum <em>F</em>(<em>e<sup>jθ</sup></em>). We can consider them as the values of <em>F</em>(<em>e<sup>jθ</sup></em>) for <em>θ </em>= 2<em>πkN</em>, since the first <em>N </em>terms in the expression for <em>F</em>(<em>e<sup>j</sup></em><sup>2<em>πk</em></sup><em><sup>/N</sup></em>) yield

<em>F</em>(<em>e</em><em>j</em>2<em>πk</em><em>/</em><em>N</em>) = +<em>…</em>

which is the same as the DFT expression. What about the remaining terms? They are repetitions of the above sum and help build up the delta function that is needed to take us from a continuous transform to discrete impulses.

<ul>

 <li>What this means is that the DFT is a sampled version of the DTFT, which is the digital version of the analog Fourier Transform</li>

</ul>

In this assignment, we want to explore how to obtain the DFT, and how to recover the analog Fourier Tranform for some known functions by the proper sampling of the function.

<h1>The DFT in Python</h1>

There are two commands in Python, one to compute the forward fourier transform and the other to compute the inverse transform. They are

numpy.fft.fft() numpy.fft.ifft()

If you import pylab, both functions are imported into the local namespace. Let us try this on a random function:

from pylab import * x=rand(100) X=fft(x) y=ifft(X) c_[x,y] print abs(x-y).max()

<ul>

 <li>When you run this code, the print statement will give you a value of about 10<sup>−15</sup>. The starting vector, <em>x</em>, and the vector got by running fft and ifft on <em>x </em>are essentially the same vector.</li>

 <li>There is another thing to note. The command c_[x,y] tells Python to make a 100×2 matrix by using <em>x </em>and <em>y </em>as the columns. This is a quick and dirty way to print both out side by side. And what this shows is that <em>x </em>is pure real while <em>y </em>is very slightly complex. This is due to the finite accuracy of the CPU so that the ifft could not exactly undo what fft did.</li>

 <li>Note that I used a 100 point vector, but it is much better to use a number that is 2<em><sup>k</sup></em>.</li>

</ul>

Let us now take a function we know about: y=sin(<em>x</em>). Let us use the fft and see what spectrum we get. We know that

<em>e</em><em>jx </em>−<em>e</em>−<em>jx</em>

<em>y </em>= sin(<em>x</em>) =

2<em>j</em>

So the expected spectrum is

1

<em>Y</em>(<em>ω</em>) =      [<em>δ </em>(<em>ω</em>−1)−<em>δ </em>(<em>ω </em>+1)]

2<em>j</em>

What do we get? We will use 128 points between 0 and 2<em>π</em>.

from pylab import * x=linspace(0,2*pi,128) y=sin(5*x) Y=fft(y) figure() subplot(2,1,1) plot(abs(Y),lw=2) grid(True) subplot(2,1,2) plot(unwrap(angle(Y)),lw=2) grid(True) show()

<ul>

 <li>h<em>eg1 </em>3i≡ from pylab import * x=linspace(0,2*pi,128) y=sin(5*x) Y=fft(y) figure() subplot(2,1,1) plot(abs(Y),lw=2) ylabel(r”$|Y|$”,size=16) title(r”Spectrum of $sin(5t)$”) grid(True) subplot(2,1,2) plot(unwrap(angle(Y)),lw=2) ylabel(r”Phase of $Y$”,size=16) xlabel(r”$k$”,size=16) grid(True) savefig(“fig9-1.png”) show()

  <ul>

   <li>We get spikes as expected. But not where we expected. There is energy at nearby frequencies as well.</li>

   <li>The spikes have a height of 64, not 0<em>.</em> We should divide by <em>N </em>to use it as a spectrum.</li>

   <li>The phase at the spikes have a phase difference of <em>π</em>, which means they are opposite signs, which is correct</li>

   <li>The actual phase at the spikes is near but not exactly correct.</li>

   <li>We haven’t yet got the frequency axis in place</li>

  </ul></li>

</ul>

What could be going wrong? The problem is that the DFT treats the position axis as another frequency axis. So it expects the vector to be on the unit circle starting at 1. Our position vector started at 0 and went to 2<em>π</em>, which is correct. The fft gave an answer in the same value. So we need to shift the <em>π </em>to 2<em>π </em>portion to the left as it represents negative frequency. This can be done with a command called <em>fftshift</em>.

The second thing that is wrong is that both 0 and 2<em>π </em>are the same point. So our 128 points need to stop at the poing just before 2<em>π</em>. The best way to do this is to create a vector of 129 values and drop the last one: x=linspace(0,2*pi,129);x=x[:-1]

The <em>ω </em>axis has a highest frequency corresponding to <em>N</em>. This is because we have <em>N </em>points between 0 and 2<em>π</em>. So ∆<em>x </em>= 2<em>π</em><em>/N</em>, and the sampling frequency becomes <em>N</em><em>/</em>2<em>π</em>. Thus <em>ω<sub>max </sub></em>= <em>N</em>. The actual frequencies go upto half that frequency only.

<ul>

 <li>h<em>eg2 </em>4i≡ from pylab import * x=linspace(0,2*pi,129);x=x[:-1] y=sin(5*x) Y=fftshift(fft(y))/128.0 w=linspace(-64,63,128) figure() subplot(2,1,1) plot(w,abs(Y),lw=2) xlim([-10,10]) ylabel(r”$|Y|$”,size=16) title(r”Spectrum of $sin(5t)$”) grid(True) subplot(2,1,2) plot(w,angle(Y),’ro’,lw=2) ii=where(abs(Y)&gt;1e-3) plot(w[ii],angle(Y[ii]),’go’,lw=2) xlim([-10,10]) ylabel(r”Phase of $Y$”,size=16) xlabel(r”$k$”,size=16) grid(True) savefig(“fig9-2.png”) show()</li>

</ul>

The plot is now as follows:

<ul>

 <li>Things have improved! The peaks are exactly at 0<em>.</em>5, and the remaining values are zero to machine precision.</li>

 <li>Better, their position along the <em>x </em>axis is correct as our function is sin(5<em>x</em>) and we expect a spike at <em>ω </em>= ±5</li>

 <li>Only the phase at the peak is meaningful and it is <em>π</em><em>/</em>2 and −<em>π</em><em>/</em>2 for the two peaks. The peak for <em>ω </em>= 5 should be 0<em>.</em>5<em>/j</em>, or 0<em>.</em>5exp(−<em>jπ</em><em>/</em>2).</li>

</ul>

Suppose we now look at AM modulation. The function we want to analyse is

<em>f</em>(<em>t</em>) = (1+0<em>.</em>1cos(<em>t</em>))cos(10<em>t</em>)

Again, we use 128 points and generate the spectrum using the above code, only changing the definition of <em>f</em>(<em>t</em>).

6 h<em>eg3 </em>6i≡ from pylab import * t=linspace(0,2*pi,129);t=t[:-1] y=(1+0.1*cos(t))*cos(10*t) Y=fftshift(fft(y))/128.0 w=linspace(-64,63,128) figure() subplot(2,1,1) plot(w,abs(Y),lw=2) xlim([-15,15]) ylabel(r”$|Y|$”,size=16)

title(r”Spectrum of $left(1+0.1cosleft(tright)right)cosleft(10tright)$”) grid(True)

subplot(2,1,2) plot(w,angle(Y),’ro’,lw=2) xlim([-15,15]) ylabel(r”Phase of $Y$”,size=16) xlabel(r”$omega$”,size=16) grid(True) savefig(“fig9-3.png”) show()

And we get …

Not what we expected!

<ul>

 <li>There are three non-zero amplitudes on both sides, as expected.</li>

 <li>All have zero phase, since we have cosines.</li>

 <li>But not three spikes. Just a broader single spike</li>

</ul>

What went wrong? Well, we did not allow for enough frequencies. Let us use 512 samples instead. But how?

<ul>

 <li>If we put more points between 0 and 2<em>π</em>,</li>

</ul>

t=linspace(0,2*pi,513);t=t[:-1]

it will give the same spacing but a high sampling frequency.

<ul>

 <li>What we need to do is to stretch the <em>t </em></li>

</ul>

t=linspace(-4*pi,4*pi,513);t=t[:-1]

This will give us tighter spacing between frequency samples, but the same time spacing – i.e., the same sampling frequency.

8 h<em>eg4 </em>8i≡ from pylab import * t=linspace(-4*pi,4*pi,513);t=t[:-1] y=(1+0.1*cos(t))*cos(10*t) Y=fftshift(fft(y))/512.0 w=linspace(-64,64,513);w=w[:-1]

figure() subplot(2,1,1) plot(w,abs(Y),lw=2) xlim([-15,15]) ylabel(r”$|Y|$”,size=16) title(r”Spectrum of $left(1+0.1cosleft(tright)right)cosleft(10tright)$”) grid(True) subplot(2,1,2) plot(w,angle(Y),’ro’,lw=2) xlim([-15,15]) ylabel(r”Phase of $Y$”,size=16) xlabel(r”$omega$”,size=16) grid(True) savefig(“fig9-4.png”) show()

Finally we are there:

Perfect! We have three clear spikes. The phases at those spikes are zero to machine precision. The location of the spikes are 9, 10 and 11 radians per sec. The height of the side bands come from

0<em>.</em>

All the amplitudes are real and positive as seen in the phase plot.

<h1>Assignment</h1>

<ol>

 <li>Work through the examples above.</li>

 <li>Generate the spectrum of sin<sup>3</sup><em>t </em>and cos<sup>3</sup><em>t</em>. Compare with what is expected.</li>

 <li>Generate the spectrum of cos(20<em>t </em>+5cos(<em>t</em>)). Plot phase points only where the magnitude is significant (say greater than 10<sup>−3</sup>). What do you think is happening?!</li>

 <li>The Gaussian exp is not “bandlimited” in frequency. We want to get its spectrum accurate to 6 digits. Try different time ranges, and see what gets you a frequency domain that is so accurate. (The Fourier Transform of a Gaussian is a</li>

</ol>

Gaussian in <em>ω</em>. Look up the transform pair to confirm)