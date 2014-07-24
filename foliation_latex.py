r"""
General supporting math functions.

AUTHORS:

- Balazs Strenner (2014-06-16): initial version

EXAMPLES::


"""

#*****************************************************************************
#       Copyright (C) 2014 Balazs Strenner <strennerb@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.functional import N

from matplotlib.colors import ColorConverter
from base import *

_tikzcolors = ["red", "green", "blue", "cyan", "magenta", "yellow", 
    "gray", "brown", "lime", "olive", "orange", "pink", "purple", 
    "teal", "violet"]

def _tikzcolor(n):
    return _tikzcolors[n % len(_tikzcolors)]

draw_separatrix_command = """% command for drawing separatrices
\\newcommand\\drawSeparatrices[1]{
\\foreach \\x/\\yone/\\ytwo in #1{
\\draw (\\x,\\yone) -- (\\x,\\ytwo);
}
}

"""

coloring_str = """% coloring

\\foreach \\list in \\endpointarray{
\\def\\lastx{}
\\def\\lastsign{}

\\foreach \\x/\\sign/\\notused/\\col/\\isflipped/\\mp[count=\\i from 0,\
remember=\\x as \\lastx, remember=\\sign as \\lastsign] in \\list{
\\ifnum \\i > 0

\\def\\darkcolor{\\col!\\colorstrength!white}
\\ifnum \\isflipped = 0 % interval is not flipped
\\def\\bp{\\colorstrength}
\\def\\ep{0}
\\else
\\def\\bp{0}
\\def\\ep{\\colorstrength}
\\fi

\\ifnum \\wrappingindex = \\i

\\ifnum \\sign = -1
\\fill[left color=\\col!\\bp!white, right color=\
\\col!\\mp!white]  (\\lastx, 0) rectangle (1,\\lastsign*0.5);
\\fill[left color=\\col!\\mp!white, right color=\
\\col!\\ep!white]  (0,0) rectangle (\\x, \\sign*0.5);
\\else
\\fill[left color=\\col!\\bp!white, right color = \\col!\\ep!white] \
(\\lastx, 0) rectangle (\\x, \\sign*0.5);
\\fi

\\else

\\fill[left color=\\col!\\bp!white, right color = \\col!\\ep!white] \
(\\lastx, 0) rectangle (\\x, \\sign*0.5);

\\fi
\\fi
}
}

"""


define_pos_str = r"""% Defines \pos as "above" or "below" whether #1 is 1 or -2
\newcommand\definepos[1]{
\ifnum #1 = 1
\def\pos{above}
\else
\def\pos{below}
\fi
}

"""

separatrices_str = """% separatrices
\\foreach \\list in \\endpointarray{
\\foreach \\x/\\sign in \\list{
\\draw (\\x,0) -- (\\x,\\sign*0.5);
}
}

"""

horizontal_sides_str = """% horizontal sides
\\foreach \\sign in {-1,1}{
\\draw[dashed] (0,\\sign*0.5) -- +(1,0);
}

"""

center_line_str = """% center line
\\draw[very thick] (0,0) -- (1,0);

"""

singularities_str = """% singularities
\\foreach \\list in \\endpointarray{
\\foreach \\x/\\sign/\\singcol in \\list{
\\filldraw[fill=\\singcol, draw=black] (\\x,\\sign*0.5) circle (\\rad);
}
}

"""

labels_str = """% labels
\\foreach \\x/\\sign/\\label in \\midpointarray{
\\definepos{\\sign}
\\node at (\\x,0) [\\pos] {\\label};
}

"""

lengths_str = """% lengths
\\foreach \\x/\\sign/\\label/\\length in \\midpointarray{
\\definepos{\\sign}
\\node at (\\x,\\sign*0.5) [\\pos] {\\length};
}

"""

train_track_str = r"""
% Arrow at a certain position of an edge, 0.5 is the center.
\tikzset{
->-/.style={
decoration={markings,
mark= at position #1 with {\arrow{>}},},
postaction={decorate}}}

% Drawing half of a center edge.
% The starting and ending points are (#1,#2) and (#3,#4).
% #5 is the label. #6 is drawing options, e.g. arrow.
\newcommand\drawedge[6]{
\pgfmathsetmacro\diff{(#4-#2)/2}
\draw[#6] (#1,#2) .. controls +(0,\diff) and +(0,-\diff) .. (#3,#4)
node[midway,\pos\space left] {#5};
}

% Drawing the train track.
% #1 - list of data for each interval
% #2 - list of codings of transformations that indentify the two
% vertical sides and the identity
\newcommand\traintrack[2]{
\foreach \x/\sign/\mid/\final/\outgoing/\ei[count=\i] in #1{
\definepos{\sign}

% the vertices v_i
\fill (\x, \sign*0.25) circle (\rad) node[right] {$v_\i$};

\ifnum \outgoing = 1
\def\starty{0.25}
\def\endy{0.5}
\def\labl{\ei}
\else
\def\starty{0.5}
\def\endy{0.25}
\def\labl{}
\fi

% the edges e_i
\draw[->-=0.5] (\x,\sign*\starty) -- node[midway,right] {\labl} (\x,\sign*\endy);

\clip (0,-0.5) rectangle (1,0.5);

% the edges f_i
\foreach \xshift/\yscale in #2{
\begin{scope}[xshift=\xshift cm,yscale=\yscale]
\drawedge{\x}{\sign*0.25}{\mid}{0}{$f_{\i}$}{->-=0.5}
\drawedge{\mid}{0}{\final}{-\sign*0.25}{}{}
\end{scope}
}
}
}

"""


class FoliationLatex(SageObject):
    """
    Options:

        - ``size`` - the width of the Tikz picture in cm's. (Usually the height
          as well, though that can be smaller.

        - ``color_strength`` - an integer on a scale of 0 to 100. The lower
          the strength the more white is mixed with the colors. 0 means no
          coloring.

        - ``interval_labelling`` - boolean, whether the intervals are 
          labelled

        - ``length_labelling`` - boolean, whether the lengths of the intervals
          are indicated on the picture. For simpler Foliations this could be
          useful, but if there are really short intervals, it is better to 
          turn this off since there is no room for these decimals.
    """

    __default_options = {
        'scale':15,
        'color_strength':50,
        'interval_labelling':True,
        'length_labelling':True,
        'separatrix_color':'darkgreen',
        'separatrix_draw_options':'dashed,thick',
        'transverse_curve_color':'red',
        'transverse_curve_draw_options':'very thick'
    }

    def __init__(self, foliation, **options):
        self._foliation = foliation
        self._options = {}
        self._separatrices = []
        self.set_options(**options)

    def _repr_(self):
        return repr(self._options)

    def set_option(self, option_name, option_value = None):
        if option_value == None:    # clear the option, if set
            if option_name in self._options:
                del self._options[option_name]
        else:
            self._options[option_name] = option_value

    def set_options(self, **kwds):
        if kwds:
            for name, value in kwds.items():
                self.set_option(name, value)

    def get_option(self, option_name):
        if not(option_name in self.__default_options):
            raise ValueError( "%s is not a Latex option for Foliation." % option_name )
        if option_name in self._options:
            return self._options[option_name]
        else:
            return self.__default_options[option_name]


    def tikz_picture(self, separatrices = [], train_tracks = [],
                     transverse_curves = []):
        r"""
        Returns the Latex/Tikz representation of the foliation.

        INPUT:

        - ``separatrices`` -- 

        OUTPUT:

        - string - the latex representation

        """
        from sage.misc.latex import latex

        scale_size = self.get_option('scale')
        color_strength = self.get_option('color_strength')
        interval_labelling = self.get_option('interval_labelling')
        length_labelling = self.get_option('length_labelling')

        latex.add_to_preamble('\\usepackage{tikz}\n')
        latex.add_to_preamble('\\usetikzlibrary{decorations.markings}\n')
        latex.add_to_preamble('\\usetikzlibrary{arrows}\n')

        fol = self._foliation
        
        s = '\\begin{{tikzpicture}}[scale = {0},>=stealth\','\
            'font=\\tiny]\n'.format(scale_size)
        s += '\\def\\rad{0.005}\n'
        s += '\\def\\colorstrength{{{0}}}\n'.format(color_strength)

        if not fol.is_bottom_side_moebius():
            wrappingindex = fol.num_intervals(BOTTOM)
        else:
            for interval in fol.intervals():
                if interval.is_wrapping(fol):
                    wrappingindex = interval.index + 1
                    break
                
        s += '\\def\\wrappingindex{{{0}}}\n\n'.format(wrappingindex)            
        s += '\\def\\midpointarray{'
        
        for interval in fol.intervals():
            signed_label = str(interval.label(fol))
            if interval.is_flipped(fol):
                signed_label = '-' + signed_label

            r = '{mid}/{side}/{label}/{length}, '
            s += r.format(mid = round(interval.endpoint(MID, fol),7),
                          side = (-1)**interval.endpoint_side(MID, fol),
                          label = signed_label,
                          length = round(interval.length(fol), 4))
            
            
        s = s[:-2] + '}\n\n\\def\\endpointarray{{0/1/'
        s += _tikzcolor(fol.intervals()[0].which_singularity(fol))

        for interval in fol.intervals():
            if interval == (BOTTOM,0):
                s += '}, {' + endpoint_entry(interval.prev(fol), fol,
                                             color_strength)
            s += ', ' + endpoint_entry(interval, fol,
                                       color_strength)


        s += '}}\n\n'
        s += coloring_str + define_pos_str + separatrices_str + \
             horizontal_sides_str + center_line_str + singularities_str

        if interval_labelling:
            s += labels_str

        if length_labelling:
            s += lengths_str

        cc = ColorConverter()

        if len(separatrices) > 0 or len(transverse_curves) > 0:
            s += draw_separatrix_command
        
        if len(separatrices) > 0:
            rgb = cc.to_rgb(self.get_option('separatrix_color'))
            s += '% separatrices\n'
            s += '\\definecolor{{separatrixcolor}}{{rgb}}{{{0},{1},{2}}}\n'.format(*rgb)
            opts = self.get_option('separatrix_draw_options')
            s += '\\tikzstyle{{separatrix opts}}=[{opts}]\n'.format(opts = opts)
            s += tikz_of_separatrices(separatrices)

        # for separatrix in separatrices:
        #     s += tikz_of_separatrix(separatrix)

        if len(transverse_curves) > 0:
            rgb = cc.to_rgb(self.get_option('transverse_curve_color'))
            s += '% transverse curves'
            s += '\\definecolor{{curvecolor}}{{rgb}}{{{0},{1},{2}}}\n'.format(*rgb)
            opts = self.get_option('transverse_curve_draw_options')
            s += '\\tikzstyle{{curve opts}}=[{opts}]\n'.format(opts = opts)
            
        for curve in transverse_curves:
            s += tikz_of_curve(curve)

        if len(train_tracks) > 0:
            s += train_track_str
            
        for train_track in train_tracks:
            s += tikz_of_train_track(train_track, fol)

        s += '\\end{tikzpicture}\n'
        return s

def endpoint_entry(interval, fol, color_strength):
    begin_percent = color_strength
    end_percent = 0
    middle_percent = ''

    if interval.is_flipped(fol) and interval > interval.pair(fol):
        begin_percent, end_percent = end_percent, begin_percent

    r = '{x}/{side}/{singcol}/{col}/{isflipped}/{mp}'
    rep = fixed_rep(interval, fol)
    side = (-1)**interval.endpoint_side(RIGHT, fol) \
           if rep < 1 else (-1)**interval.endpoint_side(LEFT,
                                                        fol)
    lep = interval.endpoint(LEFT, fol)
    if rep < lep:
        middle_percent = color_strength * rep / (1 + rep - lep)
        if begin_percent == 0:
            middle_percent = color_strength - middle_percent

    return r.format(x = rep if rep > 0 else 1,
                    side = side,
                    singcol = _tikzcolor(interval.next(fol).\
                                         which_singularity(fol)),
                    col = _tikzcolor(interval.numerical_label(fol)),
                    isflipped = 1 if interval.is_flipped(fol) else 0,
                    mp = middle_percent)

    

def fixed_rep(interval, fol):
    rep = round(interval.endpoint(RIGHT, fol),7)
    return rep if rep > 0 else 1
    

def shift_point(x, direction = None):
    if direction == None:
        return x
    SHIFT_BY = 0.005
    return x - SHIFT_BY * (-1) ** direction


def get_y(side):
    return (-1)**side * 0.5



def tikz_of_train_track(tt, fol):
    s = '\\def\\ttdata{'
    shift = 0.5 if fol.is_bottom_side_moebius() else 0.0
    fac = 1 if fol.is_bottom_side_moebius() else 0.5

    transformations = '0/1,1/-1,-1/-1' if fol.is_bottom_side_moebius() \
                      else '0/1,1/1,-1/1'
    
    for interval in fol.intervals():
        if interval.pair(fol) > interval:
            ei = '$e_{0}$'.format(interval.numerical_label(fol) + 1)
            outgoing = 1
        else:
            ei = ''
            outgoing = 0
        edge = tt.get_edge_from(interval, 'center')
        x, final = [round(edge[i].endpoint(MID, fol),7) for i in [0,1]]
        sign = (-1)**edge[0].endpoint_side(MID, fol)

        overlap_length = min(mod_one(edge[END].raw_endpoint(RIGHT, fol) + shift -
                                     edge[START].raw_endpoint(LEFT, fol)),
                             mod_one(edge[START].raw_endpoint(RIGHT, fol) -
                                     edge[START].raw_endpoint(LEFT, fol)))

        # the midpoint calculated from START
        mid = round(x - edge[START].length(fol)*fac +\
                    overlap_length*fac,7)

        # the midpoint calculated from END
        mid2 = round(final - edge[END].length(fol)*fac +
                     2*mod_one(edge[START].raw_endpoint(LEFT, fol) + shift -
                               edge[END].raw_endpoint(LEFT, fol))*fac +\
                     overlap_length*fac,7)

        # if they are different, then the edge crosses the vertical
        # boundaries, so 'final' has to be shifted
        final += mid - mid2

        s += '{x}/{sign}/{mid}/{final}/{outgoing}/{ei}'.\
             format(x = x, sign = sign, mid = mid, final = final,
             outgoing = outgoing, ei = ei) + ', '

    s = s[:-2] + '}\n\\def\\transformations{' +\
        transformations + '}\n'

    return s + '\\traintrack{\\ttdata}{\\transformations}\n\n'





def draw_segment(u, v):
    return '\\draw[color=curvecolor,curve opts] ({u},0) -- ({v},0);\n'.\
        format(u = u, v = v)

def tikz_of_curve(transverse_curve):
    arcshift = RIGHT if transverse_curve.direction() == LEFT else LEFT
    ep = [transverse_curve.arc()[i] for i in [0,1]]
    u, v = [shift_point(ep[i], arcshift) for i in [0,1]]

    if u < v:
        s = draw_segment(u, v)
    else:
        s = draw_segment(u, 1)
        if v >= 0:
            # otherwise we don't want the little segment between 0 and v
            s += draw_segment(0, v)


    sep = [transverse_curve.separatrix(i) for i in [0,1]]
    hshifts = [(arcshift + 1) % 2 if sep[i].is_flipped() else arcshift
               for i in [0,1]]
    # sep_shifts = 
    # for i in range(2):
    #     s += tikz_of_separatrix(transverse_curve.separatrix(i),
    #                                   hshifts[i])

    return tikz_of_separatrices(sep, hshifts) + s + '\n\n'
                                

def tikz_of_separatrices(separatrices, hshifts = None):
    color = 'separatrixcolor' if hshifts == None else 'curvecolor'
    opts = 'separatrix opts' if hshifts == None else 'curve opts'
    
    int_entries = []
    for i in range(len(separatrices)):
        hshift = hshifts[i] if hshifts else None
        int_entries.extend(intersection_entries(separatrices[i], hshift))
    s = make_intersection_list(int_entries)
    s += """

    \\begin{{scope}}[color={col},{opt}]
    \\drawSeparatrices{{\\intersections}}
    \\end{{scope}}

    """.format(col = color, opt = opts)
    return s

    


def intersection_entries(separatrix, hshift = None):
    entries = []
    fol = separatrix.foliation
    for i in range(separatrix.num_intersections() - 1):
        x = separatrix.get_intersection(i)
        x = shift_point(x, hshift)
        entries.append('{x}/-0.5/0.5'.format(x = x))
        if separatrix.get_tt_edge(2 * i + 1).start().is_flipped(fol):
            hshift = (hshift + 1) % 2 if hshift != None else None

    end_x = shift_point(round(separatrix.endpoint,7), hshift)
    end_y = get_y(separatrix.end_side())
    entries.append('{x}/0/{y}'.format(x=end_x, y = end_y))

    return entries

def make_intersection_list(intersection_entry_list):
    return "\\def\\intersections{{{0}}}".\
        format(', '.join(intersection_entry_list))


