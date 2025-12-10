import numpy as np
import scipy.integrate as integrate
from scipy.special import hermite
import math
from matplotlib import pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import random as rand
import itertools

# ===== WAVEFORM/CURVE CONSTANTS =====

# Quantum system parameters
HILBERT_DIM = 4  # Increased to 4 basis states for more complexity

# Wave function parameters
WAVE_COEFFICIENTS = [0.4+0.2j, 0.3-0.4j, 0.2+0.3j, 0.1-0.1j]  # Complex coefficients for richer interference patterns

# Animation parameters
ANIMATION_FRAMES = 3600  # Total number of animation frames
TIME_FACTOR = 0.8        # Time evolution speed multiplier (slower for complex patterns)
SAMPLING_POINTS = 120    # Number of spatial sampling points (higher resolution)

# Physical potential
POTENTIAL_FUNCTION = lambda x: 1/2*x**2  # Harmonic oscillator: V(x) = (1/2)x²

# Visualization limits
AXIS_X_LIMITS = (-6, 6)    # Position axis range (expanded for more complex patterns)
AXIS_Y_LIMITS = (-2, 2)    # Imaginary part axis range (expanded for higher amplitudes)
AXIS_Z_LIMITS = (-2, 2)    # Real part axis range (expanded for higher amplitudes)

# Trail effect parameters
TRAIL_COUNT = 16           # Number of trail copies
TRAIL_OFFSET_STEP = 32     # Frame offset between trails
TRAIL_FADE_BASE = TRAIL_COUNT + 1  # Base for alpha calculation (ensures proper fading)

# Color parameters for trail fading
TRAIL_BACKGROUND_COLOR = 30  # RGB value for faded trail background (0-255)
ORIGINAL_BLUE_VALUE = 255    # Original blue channel value
ORIGINAL_GREEN_VALUE = 255   # Original green channel value

# Camera parameters
CAMERA_ELEVATION = 30     # Camera elevation angle (degrees)
CAMERA_AZIMUTH = -60      # Camera azimuth angle (degrees)

# Animation timing
ANIMATION_INTERVAL = 10   # Milliseconds between frames
TIME_SCALING_FACTOR = 100  # Divisor for time scaling in animation

# Simple numerical derivative function to replace scipy.misc.derivative
def derivative(func, x, dx=1e-6, n=1):
    """Numerical derivative using finite differences"""
    if n == 1:
        return (func(x + dx) - func(x - dx)) / (2 * dx)
    elif n == 2:
        return (func(x + dx) - 2*func(x) + func(x - dx)) / (dx**2)
    else:
        raise ValueError("Only n=1 and n=2 supported")

# Quantum mechanics simulator for (initially) finite dimensional
# approximations of Hilbert spaces provided a basis.

'''
Code plan and structure:

One-dimensional:
1. Given an arbitrary Hamiltonian H, find the eigenvalues E_n.
2. Estimate the matrix form of the Hamiltonian in a finite-dimensional 
   approximation in the QHO basis.
  a. Determine the numerical solutions to the differential equations H.psi = E_n.psi.
     (Call it numSolution[n])
  b. Perform functional optimization over the parameter space by finding the set
     of coefficients, [c_0, c_1, ..., c_n] such that the value
     
     QHOBasisApprox = numSolution[n](x) - sum([c[n] * hilbert.eigenbasis(n, x)
     for n in range(hilbert.dim)])

     assumes its minimum value.
  c. Take the basis approximation and form an approximation of the Hamiltonian,

     H = np.empty(2*[len(numSolution])
     for i, j in itertools.product(range(hilbert.dim), range(hilbert.dim)):
         H[i][j] = lambda x: np.conj(QHOBasisApprox[i](x)) * QHOBasisApprox[j](x)
  d. Set the unitary time operator to be

     U = lambda x, t: np.exp(-i/hbar * H(x) * t)
  e. Evolve the wave function to be
     psi = U(x, t) * initWaveFunc


Three-dimensional:
'''

class FunctionSampler:
    # Secret occult tricks to make computing our basis functions cheap.

    def __init__(
        self,
        f,
        minX,
        maxX,
        numSamples,
    ):
        self.minX = minX
        self.maxX = maxX
        self.range = maxX - minX
        self.numSamples = numSamples
        
        self.domain = np.linspace(
            minX,
            maxX,
            numSamples,
        )
        self.image = [
            f(x)
            for x in self.domain
        ]

    def __call__(self, x):
        x = min(
            max(
                x,
                self.minX,
            ),
            self.maxX - 1
        )
        x -= self.minX
        x = round(
            x * (self.numSamples / self.range)
        )

        return self.image[int(x)]

class HilbertSpace:

    def __init__(
        self,
        dim,
        basis = 'QHO',
        hamiltonianPotential = lambda x: x**4,
        approxMethod = 'QHO',
    ):
        self.dim = dim

        if hamiltonianPotential != None:
            self.V = hamiltonianPotential

        if basis == 'QHO' or approxMethod == 'QHO':
            self.QHOBasis = lambda n, x: (
                1/np.sqrt(
                    2**n * math.factorial(n)
                ) * np.pi**(1/4)) * np.exp(-x**2/2) * hermite(n)(x)
            QHOBasisApprox = []
            for n in range(self.dim):
                QHOBasisApprox.append(
                    FunctionSampler(
                        lambda x: self.QHOBasis(n,x),
                        -15,
                        15,
                        2000,
                    )
                )

            self.QHOBasis = lambda n, x: QHOBasisApprox[n](x)

        elif basis == 'Square Well' or approxMethod == 'Square Well':
            self.SWBasis = lambda n, x: (
                np.sqrt(2)
                * np.sin(
                    (n + 1)
                    * np.pi
                    * (x + 1 / 2)
                )
                if abs(x) < 1 / 2
                else 0
            )

        if approxMethod == 'QHO':
            self.H = np.empty(
                [
                    self.dim,
                    self.dim,
                ]
            )
            for i, j in itertools.product(
                range(self.dim),
                range(self.dim),
            ):
                HFunc = lambda x: np.conj(
                    self.QHOBasis(i, x)
                ) * (
                    -1/2 * derivative(
                        lambda xPrime: self.QHOBasis(
                            j,
                            xPrime,
                        ),
                        x,
                        dx=1e-6,
                        n=2
                    )
                    + self.V(x) * self.QHOBasis(j,x)
                )
                self.H[i][j] = integrate.quad(HFunc, -np.inf, np.inf)[0]
                
        if basis == 'QHO':
            self.eigenbasis = self.QHOBasis
            self.eigenvalues = lambda n: 1/2 + n
            
        elif basis == 'Square Well':
            self.eigenbasis = self.SWBasis
            self.eigenvalues = lambda n: (n+1)**2 * np.pi**2
            
class WaveFunction:

    def __init__(
        self,
        hilbertSpace,
        initWaveFunc = None,
        coeff = None,
    ):
        self.hilbert = hilbertSpace

        if coeff != None:
            self.coeff = coeff
                
        else:
            self.evaluate = lambda x: initWaveFunc(x)
            self.coeff = self.orthogonalBasisProjection(initWaveFunc)

        evaluate = lambda x, t: sum(
            [
                self.coeff[n] * self.hilbert.eigenbasis(n, x) * self.phaseFactor(n, t)
                for n in range(self.hilbert.dim)
            ]
        )
        normF = self.normalize(lambda x: evaluate(x, 0))
        
        self.evaluate = lambda x, t: evaluate(x, t) / normF
        
    def phaseFactor(self, n, t):
        return np.exp(-1j * self.hilbert.eigenvalues(n) * t)

    def normalize(self, waveFunc):
        return np.sqrt(integrate.quad(
            lambda x: np.absolute(
                waveFunc(x)
            )**2,
            -np.inf,
            np.inf,
        )[0])

    def orthogonalBasisProjection(self, waveFunc):
        coeff = []
        for n in range(hs.dim):
            coeff.append(
                integrate.quad(
                    lambda x: np.conj(
                        self.hilbert.eigenbasis(n, x)
                    ) * waveFunc(x),
                    -np.inf,
                    np.inf,
                )[0])

        return coeff

class Plotter:

    def __init__(self):
        pass

    def plotWaveFunction2d(
        self,
        waveFunction,
        samples = 50,
        frames = 3600,
        timeFactor = 1.0/32.0,
        saveName = None,
    ):

        fig = plt.figure()
        ax1 = plt.axes(xlim=(-5, 5), ylim=(-1,1))
        line, = ax1.plot([], [], lw=2)
        plt.xlabel('$x$')
        plt.ylabel(r'$\psi(x)$')
        plt.xticks([],[])
        plt.yticks([],[])

        plotlays, plotcols = [2], ["blue","red","green"]
        lines = []
        for index in range(3):
            lobj = ax1.plot(
                [],
                [],
                lw=2,
                color=plotcols[index]
            )[0]
            lines.append(lobj)

        def init():
            for line in lines:
                line.set_data([], [])
            return lines

        def animate(i):
            posValues = np.linspace(-5, 5,samples)
            amplitudeValues = [
                waveFunction.evaluate(
                    x,
                    timeFactor * i/50
                )
                for x in posValues
            ]

            reValues = np.real(amplitudeValues)
            imValues = np.imag(amplitudeValues)
            probValues = np.absolute(amplitudeValues)**2

            ylist = [reValues, imValues, probValues]

            for lnum, line in enumerate(lines):
                line.set_data(posValues, ylist[lnum])
                
            return lines

        anim = animation.FuncAnimation(
            fig,
            animate,
            init_func = init,
            frames = frames,
            interval = 10,
            blit = True,
        )

        if saveName != None:
            anim.save(
                saveName + '.mp4',
                fps = 30,
                extra_args = ['-vcodec', 'libx264'],
            )

        plt.show()

    def plotWaveFunction3d(
        self,
        waveFunction,
        samples=SAMPLING_POINTS,
        frames=ANIMATION_FRAMES,
        timeFactor=TIME_FACTOR,
        saveName=None,
    ):

        fig = plt.figure()
        ax1 = fig.add_axes([0,0,1,1], projection='3d', proj_type='persp')
        ax1.set_xlabel('$x$')
        ax1.set_ylabel(r'Im$(\psi)$')
        ax1.set_zlabel(r'Re$(\psi)$')
        
        ax1.set_xlim(AXIS_X_LIMITS)
        ax1.set_ylim(AXIS_Y_LIMITS)
        ax1.set_zlim(AXIS_Z_LIMITS)

        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_zticks([])

        plotlays, plotcols = [2], ["blue","red","green"]
        
        lines = []
        pts = []
        
        # Create main curves (current state) - only blue curve
        for index in [0]:  # Removed green curve (index 2)
            lobj = ax1.plot([],[], [],lw=2,color=plotcols[index])[0]
            pobj = ax1.plot([],[], [],lw=2,color=plotcols[index])[0]
            lines.append(lobj)
            pts.append(pobj)
        
        # Create trail curves (past states) - only for blue curve
        trail_lines = []
        for trail_idx in range(TRAIL_COUNT):
            alpha = 1.0 - (trail_idx + 1) / TRAIL_FADE_BASE
            for index in [0]:  # Only blue curve trails
                # Interpolate color from original to background
                if index == 0:  # blue to dark
                    r = int(0 * (1-alpha) + TRAIL_BACKGROUND_COLOR * alpha)
                    g = int(0 * (1-alpha) + TRAIL_BACKGROUND_COLOR * alpha) 
                    b = int(ORIGINAL_BLUE_VALUE * (1-alpha) + TRAIL_BACKGROUND_COLOR * alpha)
                # Removed green curve trail logic
                
                color = f'#{r:02x}{g:02x}{b:02x}'
                lobj = ax1.plot([],[], [],lw=1.5,color=color,alpha=alpha*0.8)[0]
                trail_lines.append(lobj)

        def init():
            for line, pt in zip(lines, pts):
                line.set_data([], [])
                line.set_3d_properties([])
            
            for trail_line in trail_lines:
                trail_line.set_data([], [])
                trail_line.set_3d_properties([])

            return lines + pts + trail_lines

        # Set initial camera position (disabled frame-based camera movement)
        ax1.view_init(elev=CAMERA_ELEVATION, azim=CAMERA_AZIMUTH)

        def frame_path(i):
            # Disabled: camera position no longer changes with frame
            return CAMERA_ELEVATION, CAMERA_AZIMUTH

        def animate(i):
            # Removed: elev, azim = frame_path(i)
            # Removed: ax1.view_init(elev, azim)
            
            posValues = np.linspace(
                -5,
                5,
                samples
            )
            
            # Calculate current state
            amplitudeValues = [
                waveFunction.evaluate(
                    x,
                    timeFactor * i / TIME_SCALING_FACTOR
                )
                for x in posValues
            ]

            reValues = np.real(amplitudeValues)
            imValues = np.imag(amplitudeValues)
            probValues = np.absolute(
                amplitudeValues
            ) ** 2

            ylist = [
                (imValues, reValues),  # Only blue curve: x vs Im(ψ) vs Re(ψ)
            ]

            # Update main curves
            for n, line in enumerate(lines):
                line.set_data(
                    posValues,
                    ylist[n][0]
                )
                line.set_3d_properties(
                    ylist[n][1]
                )

            # Update trail curves - only for blue curve
            trail_idx = 0
            for offset in range(TRAIL_OFFSET_STEP, TRAIL_OFFSET_STEP * TRAIL_COUNT + 1, TRAIL_OFFSET_STEP):
                past_i = max(0, i - offset)
                
                # Calculate past state
                past_amplitudeValues = [
                    waveFunction.evaluate(
                        x,
                        timeFactor * past_i / TIME_SCALING_FACTOR
                    )
                    for x in posValues
                ]

                past_reValues = np.real(past_amplitudeValues)
                past_imValues = np.imag(past_amplitudeValues)

                past_ylist = [
                    (past_imValues, past_reValues),  # Only blue curve trails
                ]

                # Update trail lines for this time offset
                for n in range(1):  # Only 1 curve type now
                    trail_lines[trail_idx].set_data(
                        posValues,
                        past_ylist[n][0]
                    )
                    trail_lines[trail_idx].set_3d_properties(
                        past_ylist[n][1]
                    )
                    trail_idx += 1

            fig.canvas.draw()
            return lines + trail_lines

        anim = animation.FuncAnimation(
            fig,
            animate,
            init_func=init,
            frames=frames,
            interval=ANIMATION_INTERVAL,
            blit=True
        )

        if saveName != None:
            anim.save(
                './video/animations/' + saveName + '.mp4',
                fps = 30,
                writer = 'avconv',
                codec = 'libx264'
            )

        plt.show()

p = Plotter()
hs = HilbertSpace(
    dim=HILBERT_DIM,
    hamiltonianPotential=POTENTIAL_FUNCTION,
    basis='QHO',
    approxMethod=None,
)
psi = WaveFunction(hs, coeff=WAVE_COEFFICIENTS)
#psi = WaveFunction(hs, lambda x: 1 if abs(x) < 5 else 0)
print('c_n = ' + str(psi.coeff))

p.plotWaveFunction3d(
    psi,
    samples=SAMPLING_POINTS,
    frames=ANIMATION_FRAMES,
    timeFactor=TIME_FACTOR,
    saveName=None,
)

'''
Ideas:

- Take a state vetor and evolve it using the unitary evolution operator defined
  as the exponentiated Hamiltonian. Recall that the Hamiltonian defines the 
  energy eigenbasis. If you have an eigenbasis and eigenvalues, simply define
  the Hamiltonian according to the spectral theorem. If you instead have a
  Hamiltonian first, we have to decompose it into a finite-dimensional
  approximation.

  One method to do this may be to calculate the first n eigenstates and
  eigenvalues of the Hamiltonian, and then to use spectral theorem to
  simply define the Hamiltonian.

- Long term goals may include determining the unitary evolution of a state in a
  time-varying potential, determining the projection onto the eigenstate of an
  observable O (to allow for measurements, which can be determined by sampling
  our probability distribution in the corresponding observable picture)

- Mapping the position wave function to the corresponding momentum picture using
  Fourier transforms, and other methods?

- Checking Vachaspati's hypothesis that the expectation value of a wave function is
  related to the RMS value of a classical wave??? Not sure what he was thinking,
  perhaps ask JJ.

- Coherent states.
'''
