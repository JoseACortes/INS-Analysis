def make_window(wave, window_min, window_max):
    """
    Make a window from a wave and a window dict.
    """
    window_wave = wave[:, 
        (wave[0] >= window_min) & (wave[0] <= window_max)
        ]
    return window_wave