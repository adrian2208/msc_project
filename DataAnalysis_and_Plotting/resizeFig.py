#ALL CREDIT GOES TO https://jwalton.info/Embed-Publication-Matplotlib-Latex/
def set_size(width = 373.45274, fraction=1):
    fig_width_pt = width * fraction
    inches_per_pt = 1 / 72.27
    golden_ratio = (5**.5 - 1) / 2
    fig_width_in = fig_width_pt * inches_per_pt
    fig_height_in = fig_width_in * golden_ratio
    fig_dim = (fig_width_in, fig_height_in)
    return fig_dim
