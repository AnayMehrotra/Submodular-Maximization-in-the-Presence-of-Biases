home_folder = "."

rcParams.update({
        'text.usetex': False,
        'font.family': 'stixgeneral',
        'mathtext.fontset': 'stix',
        'figure.figsize': (10,7),
})

def file_str():
    """ Auto-generates file name."""
    now = datetime.datetime.now()
    return now.strftime("H%HM%MS%S_%m-%d-%y")

rand_string = lambda length: ''.join(random.choice(string.ascii_lowercase) for i in range(length))

def pdf_savefig():
    """ Saves figures as pdf """
    tmptmp = file_str()+rand_string(5)
    print(home_folder+f"/figs/{tmptmp}.pdf")
    plt.savefig(home_folder+f"/figs/{tmptmp}.pdf")

def eps_savefig():
    """ Saves figure as encapsulated postscript file (vector format)
        so that it isn't pixelated when we put it into a pdf. """
    pdf_savefig()
