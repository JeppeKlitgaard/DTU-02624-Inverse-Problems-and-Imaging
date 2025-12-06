def forward_heat_operator(x: npt.NDArray, T: float) -> npt.NDArray:
    # An efficient way to solve this is to turn the convolution
    # into elementwise multiplication in the frequency domain
    # In Homework 2 we used DST, but here we are specifically asked to use FFT
    # So we follow odd-expansion as in Kim's given code (ModelProblems)

    N = x.shape[0]
    assert N == x.shape[1], "Input must be square"

    # Spatial scale
    ds = np.pi / (N + 1)

    # Wavenumbers
    N_ext = 2 * (N + 1)
    frequencies = fft.fftfreq(N_ext, d=ds)
    wavenumbers = 2 * np.pi * frequencies
    n = wavenumbers
    m = wavenumbers
    n_grid, m_grid = np.meshgrid(n, m, indexing="ij")

    # Decay factors
    k2 = n_grid**2 + m_grid**2  # Combined wavenumber
    decay_factors = np.exp(-(k2) * T)  # Here the same as eigvals or singvals

    # Construct odd expansion (Pad with zeros, reflect and negate)
    x_ext = np.zeros((N_ext, N_ext))  # Pad
    x_ext[1:N+1, 1:N+1] = x  # Insert original
    x_ext[N+2:, 1:N+1] = -np.flipud(x)  # Reflect along horizontal
    x_ext[:, N+2:] = -np.fliplr(x_ext[:, 1:N+1])  # Reflect along vertical

    # Carry out convolution in frequency domain
    x_ext_hat = fft.fft2(x_ext)
    x_ext_hat_T = x_ext_hat * decay_factors
    y_ext = fft.ifft2(x_ext_hat_T)

    y = y_ext[1:N+1, 1:N+1]

    # Optionally, if we care about being normalised:
    # Scaling factor due to basis functions sin(...) not being normalised
    # scaling_factor = (np.pi / 2)**2
    # y = y * scaling_factor
    return np.real(y)
