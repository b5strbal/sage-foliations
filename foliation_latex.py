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


define_pos_str = r"""\newcommand\definepos[1]{
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

    def _draw_segment(self, u, v):
        return '\\draw[color=curvecolor,curve opts] ({u},0) -- ({v},0);\n'.\
            format(u = u, v = v)

    def _tikz_of_curve(self, transverse_curve):
        arcshift = RIGHT if transverse_curve.direction() == LEFT else LEFT
        ep = [transverse_curve.arc()[i] for i in [0,1]]
        u, v = [shift_point(ep[i], arcshift) for i in [0,1]]
        
        if u < v:
            s = self._draw_segment(u, v)
        else:
            s = self._draw_segment(u, 1)
            if v >= 0:
                # otherwise we don't want the little segment between 0 and v
                s += self._draw_segment(0, v)


        sep = [transverse_curve.separatrix(i) for i in [0,1]]
        hshifts = [(arcshift + 1) % 2 if sep[i].is_flipped() else arcshift
                   for i in [0,1]]
        for i in range(2):
            s += self._tikz_of_separatrix(transverse_curve.separatrix(i),
                                          hshifts[i])
        return s

    def _tikz_of_separatrix(self, separatrix, hshift = None):
        color = 'separatrixcolor' if hshift == None else 'curvecolor'
        opts = 'separatrix opts' if hshift == None else 'curve opts'

        fol = self._foliation
        s = ''
        for i in range(separatrix.num_intersections() - 1):
            x = separatrix.get_intersection(i)
            x = shift_point(x, hshift)
            s += '\\draw[color={col}, {opt}] ({x},-0.5) --'\
                 ' ({x},0.5);\n'.format(x=x, opt = opts,
                                        col = color)
            if separatrix.get_tt_edge(2 * i + 1).start().is_flipped(fol):
                hshift = (hshift + 1) % 2 if hshift != None else None

        # # the following correction is necessary to handle the difference
        # # between a separatrix that stops just before the Moebius band
        # # and another that stops right after it.
        # if self._foliation.is_bottom_side_moebius() and \
        #    separatrix.num_intersections() % 2 == 0:
        #     end_index = -2
        # else:
        #     end_index = -1
        end_x = shift_point(N(separatrix.endpoint), hshift)
        end_y = get_y(separatrix.end_side())
        s += '\\draw[color={col}, {opt}] ({x},0) -- ({x},{y});\n'.format(x=end_x,
                                                                         y=end_y,
                                                                         opt = opts,
                                                                         col = color)
        return s



    def _tikz_of_train_track(self, train_track):
        s = ''
        # for vertex in train_track.vertices():
        #     s += self._tikz_of_train_track_vertex(vertex)
        for edge in train_track.edges():
            s += self._tikz_of_train_track_edge(edge)
        return s

    # def _tikz_of_train_track_vertex(self, vertex):
    #     interval = vertex
    #     return '\\fill ({x},{y}) circle (0.005);\n'.format(
    #         x = interval.midpoint(), y = self._end_y(interval.side)/2)



    def _tikz_of_train_track_edge(self, edge):
        fol = self._foliation
        s = ''

            # r = '{startx}/{sign}/{midx}/{endx}'


        m = [N(edge[i].endpoint(MID, fol)) for i in [0,1]]
        y = [get_y(edge[i].endpoint_side(MID, fol))/2 for i in [0,1]]
        arrows = ['>','<']
        if edge[2] == 'pair':
            return "".join(['\\fill ({x},{y1}) circle (0.005);\n'\
                '\\draw[-{arr}-={pos}] ({x},{y1}) -- ({x},{y2});\n'.format(
                    pos = (-1)**i * 0.5, x = m[i],
                    y1 = y[i], y2 =2 * y[i],
                    arr = arrows[i]) for i in [0,1]])
                            
        clip = edge[START].is_wrapping(fol) or edge[END].is_wrapping(fol)
        if clip:
            s += '\\begin{scope}\n'\
                 '\\clip (0,-0.5) rectangle (1,0.5);\n'
        
        shift = 0.5 if fol.is_bottom_side_moebius() else 0.0

        overlap_length = min(mod_one(edge[END].raw_endpoint(RIGHT, fol) + shift -
                                     edge[START].raw_endpoint(LEFT, fol)),
                             mod_one(edge[START].raw_endpoint(RIGHT, fol) -
                                     edge[START].raw_endpoint(LEFT, fol)))

        x = []
        fac = 1 if fol.is_bottom_side_moebius() else 0.5
        x.append(N(m[0] - edge[START].length(fol)*fac + overlap_length*fac))
        x.append(N(m[1] - edge[END].length(fol)*fac +
                 2*mod_one(edge[START].raw_endpoint(LEFT, fol) + shift -
                           edge[END].raw_endpoint(LEFT, fol))*fac + overlap_length*fac))
                             
            
        transformations = [{''}, {''}]        
        for i in range(2):
            if x[i] < 0:
                if fol.is_bottom_side_moebius():
                    transformations[i].add('xshift=1cm,yscale=-1')                     
                else:
                    transformations[i].add('xshift=1cm')
            elif x[i] > 1:
                if fol.is_bottom_side_moebius():
                    transformations[i].add('xshift=-1cm,yscale=-1')                     
                else:
                    transformations[i].add('xshift=-1cm')

        # print m, x, y, overlap_length, right_endpoint
        
        s += center_edge_piece(m[0], y[0], x[0], 0, transformations[0], True)
        s += center_edge_piece(x[1], 0, m[1], y[1], transformations[1])

        if clip:
            s += '\\end{scope}\n'

        return s



        





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
        
        s = ''
        if len(train_tracks) > 0:
            for arrow in ['>', '<']:
                s += '\\tikzset{-' + arrow + '-/.style={\n'\
                     'decoration={markings,\n'\
                     'mark= at position #1 with {\\arrow{' + arrow + '}},'\
                     '},\n'\
                     'postaction={decorate}}}\n'

        s += '\\begin{{tikzpicture}}[scale = {0},>=stealth\','\
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
            s += r.format(mid = N(interval.endpoint(MID, fol)),
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
            
        if len(separatrices) > 0:
            rgb = cc.to_rgb(self.get_option('separatrix_color'))
            s += '\\definecolor{{separatrixcolor}}{{rgb}}{{{0},{1},{2}}}\n'.format(*rgb)
            opts = self.get_option('separatrix_draw_options')
            s += '\\tikzstyle{{separatrix opts}}=[{opts}]\n'.format(opts = opts)

        for separatrix in separatrices:
            s += self._tikz_of_separatrix(separatrix)

        if len(transverse_curves) > 0:
            rgb = cc.to_rgb(self.get_option('transverse_curve_color'))
            s += '\\definecolor{{curvecolor}}{{rgb}}{{{0},{1},{2}}}\n'.format(*rgb)
            opts = self.get_option('transverse_curve_draw_options')
            s += '\\tikzstyle{{curve opts}}=[{opts}]\n'.format(opts = opts)
            
        for curve in transverse_curves:
            s += self._tikz_of_curve(curve)

        for train_track in train_tracks:
            s += self._tikz_of_train_track(train_track)

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
    rep = N(interval.endpoint(RIGHT, fol))
    return rep if rep > 0 else 1
    

def shift_point(x, direction = None):
    if direction == None:
        return x
    SHIFT_BY = 0.005
    return x - SHIFT_BY * (-1) ** direction


def get_y(side):
    return (-1)**side * 0.5


def center_edge_piece(x1, y1, x2, y2, transformations = '',
                       has_arrow = False):
    r"""

    INPUT:

    - ``x1`` -- 

    - ``y1`` -- 

    - ``x2`` -- 

    - ``y2`` -- 

    """
    s = ''
    for transf in transformations:
        s += '\\draw[{arrow},{transf}] '\
             '({x1},{y1}) .. controls +(0,{tan}) and'\
             ' +(0,-{tan}) .. ({x2},{y2});\n'.format(
                 x1 = x1, x2 = x2, y1 = y1, y2 = y2,
                 transf = transf , tan = (y2 - y1)/2,
                 arrow = '->' if has_arrow else '')
    return s

