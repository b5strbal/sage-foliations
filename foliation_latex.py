
from sage.structure.sage_object import SageObject

from matplotlib.colors import ColorConverter
from foliation import mod_one
from constants import *

_tikzcolors = ["red", "green", "blue", "cyan", "magenta", "yellow", 
    "gray", "brown", "lime", "olive", "orange", "pink", "purple", 
    "teal", "violet"]

def _tikzcolor(n):
    return _tikzcolors[n % len(_tikzcolors)]




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
        arcshift = RIGHT if transverse_curve.direction() == 'left' else LEFT
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

        s = ''
        for i in range(separatrix.num_intersections() - 1):
            x = separatrix.get_intersection(i)
            x = shift_point(x, hshift)
            s += '\\draw[color={col}, {opt}] ({x},-0.5) --'\
                 ' ({x},0.5);\n'.format(x=x, opt = opts,
                                        col = color)
            if separatrix.get_tt_edge(2 * i + 1).start().is_flipped():
                hshift = (hshift + 1) % 2 if hshift != None else None

        # # the following correction is necessary to handle the difference
        # # between a separatrix that stops just before the Moebius band
        # # and another that stops right after it.
        # if self._foliation.is_bottom_side_moebius() and \
        #    separatrix.num_intersections() % 2 == 0:
        #     end_index = -2
        # else:
        #     end_index = -1
        end_x = shift_point(separatrix.endpoint, hshift)
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
        s = ''
        m = [edge[i].endpoint(MID) for i in [0,1]]
        y = [get_y(edge[i].endpoint_side(MID))/2 for i in [0,1]]
        if edge[2] == 'pair':
            return "".join(['\\fill ({x},{y1}) circle (0.005);\n'\
                '\\draw[->-={pos}] ({x},{y1}) -- ({x},{y2});\n'.format(
                pos = (-1)**i * 0.5, x = m[i], y1 = y[i], y2 =2*y[i]) for
                            i in [0,1]])
                            
        clip = edge[0].is_wrapping() or edge[1].is_wrapping()
        if clip:
            s += '\\begin{scope}\n'\
                 '\\clip (0,-0.5) rectangle (1,0.5);\n'
        
        shift = 0.5 if self._foliation.is_bottom_side_moebius() else 0.0

        overlap_length = min(mod_one(edge[1].raw_endpoint(1) + shift -
                                     edge[0].raw_endpoint(0)),
                             mod_one(edge[0].raw_endpoint(1) -
                                     edge[0].raw_endpoint(0)))
        # right_endpoint = mod_one(edge[0].endpoint(0) + overlap_length)

        x = []
        fac = 1 if self._foliation.is_bottom_side_moebius() else 0.5
        x.append(m[0] - edge[0].length()*fac + overlap_length*fac)
        x.append(m[1] - edge[1].length()*fac +
                 2*mod_one(edge[0].raw_endpoint(0) + shift -
                         edge[1].raw_endpoint(0))*fac + overlap_length*fac)
                             
            
        transformations = [{''}, {''}]        
        for i in range(2):
            if x[i] < 0:
                if self._foliation.is_bottom_side_moebius():
                    transformations[i].add('xshift=1cm,yscale=-1')                     
                else:
                    transformations[i].add('xshift=1cm')
            elif x[i] > 1:
                if self._foliation.is_bottom_side_moebius():
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

        s = ''
        if len(train_tracks) > 0:
            s += '\\tikzset{->-/.style={\n'\
                 'decoration={markings,\n'\
                 'mark= at position #1 with {\\arrow{>}},'\
                 '},\n'\
                 'postaction={decorate}}}\n'

        s += '\\begin{{tikzpicture}}[scale = {0},>=stealth\','\
            'font=\\tiny]\n'.format(scale_size)


        singularities = ''
        lines = ''
        fillings = ''
        labels = ''

        fol = self._foliation
        for interval in self._foliation.intervals():

            begin_percent = color_strength
            end_percent = 0
            middle_percent = None

            signed_label = str(interval.label())
            if interval.is_flipped():
                signed_label = '-' + signed_label
                if interval > interval.pair():
                    begin_percent, end_percent = end_percent, begin_percent

            x = [interval.endpoint(i) for i in [LEFT, RIGHT]]
            midx = interval.endpoint(MID)
            y = [get_y(interval.endpoint_side(i)) for i in [LEFT, RIGHT]]
            color = _tikzcolor(self._foliation.index_of_label(
                interval.label()))
            if x[1] == 0:
                x[1] = 1
                if fol.is_bottom_side_moebius():
                    y[1] = -0.5


            temp_lines = '\\draw ({x0},0) -- ({x0},{y0});\n'
            if x[0] < x[1]:
                temp_lines += '\\draw[dashed] ({x0},{y0}) -- ({x1},{y1});\n'
                temp_fillings = '\\shade[left color = {col}!{bp}!white, '\
                        'right color = {col}!{ep}!white] ({x0},0) rectangle '\
                        '({x1},{y1});\n'
            else:
                temp_lines += '\\draw[dashed] ({x0},{y0}) -- (1,{y0});\n'
                temp_lines += '\\draw[dashed] (0,-0.5) -- ({x1},-0.5);\n'
                middle_percent = color_strength * (x[1]) / (1 + x[1] - x[0])
                if begin_percent == 0:
                    middle_percent = color_strength - middle_percent

                temp_fillings = '\\shade[left color = {col}!{bp}!white, '\
                        'right color = {col}!{mp}!white] ({x0},0) rectangle '\
                        '(1,{y0});\n'
                temp_fillings += '\\shade[left color = {col}!{mp}!white, '\
                        'right color = {col}!{ep}!white] (0,0) rectangle '\
                        '({x1},{y1});\n'

            lines += temp_lines.format(x0=x[0],x1=x[1],y0=y[0],y1=y[1])
            fillings += temp_fillings.format(col=color,
                                             x0=x[0],x1=x[1],y0=y[0],y1=y[1],
                                             bp=begin_percent,
                                             ep=end_percent,
                                             mp=middle_percent)

            sing_color = _tikzcolor(interval.which_singularity())
            singularities += '\\filldraw[fill={col}, draw = black] ({x0},{y0}) '\
                    'circle (0.005);\n'.format(x0=x[0],y0=y[0], col = sing_color)
            midy = get_y(interval.endpoint_side(MID))
            above_or_below = 'above' if midy == 0.5 else 'below'

            if length_labelling:
                labels += '\\node at ({m},0) [{pos}] {{{label}}};\n'.\
                        format(m=midx, pos=above_or_below, label=signed_label)
            if interval_labelling:
                labels += '\\node at ({m},{y}) [{pos}] {{{length}}};\n'.\
                        format(m=midx, y=midy, pos=above_or_below, length=round(interval.length(), 4))

        lines += '\\draw (1,0) -- (1,{y});\n'\
                '\\draw[very thick] (0,0) -- (1,0);\n'.format(y=-0.5 if fol.is_bottom_side_moebius()
                                                              else 0.5)


        cc = ColorConverter()
        rgb = cc.to_rgb(self.get_option('separatrix_color'))
        s += '\\definecolor{{separatrixcolor}}{{rgb}}{{{0},{1},{2}}}\n'.format(*rgb)
        rgb = cc.to_rgb(self.get_option('transverse_curve_color'))
        s += '\\definecolor{{curvecolor}}{{rgb}}{{{0},{1},{2}}}\n'.format(*rgb)


        opts = self.get_option('separatrix_draw_options')
        s += '\\tikzstyle{{separatrix opts}}=[{opts}]\n'.format(opts = opts)

        opts = self.get_option('transverse_curve_draw_options')
        s += '\\tikzstyle{{curve opts}}=[{opts}]\n'.format(opts = opts)


        s += fillings + lines + singularities + labels
        for separatrix in separatrices:
            s += self._tikz_of_separatrix(separatrix)

        for curve in transverse_curves:
            s += self._tikz_of_curve(curve)

        for train_track in train_tracks:
            s += self._tikz_of_train_track(train_track)

        s += '\\end{tikzpicture}\n'
        return s

    # def _adjust_point(self, x):
    #     if not self._foliation.is_bottom_side_moebius():
    #         return x
    #     if x > 0.5 - epsilon:
    #         return 2 * (x - 0.5)
    #     return 2 * x


    # def _get_side(self, side, point):
    #     if not self._foliation.is_bottom_side_moebius():
    #         return side
    #     if point < 0.5:
    #         return TOP
    #     else:
    #         return BOTTOM

    # def _get_y(self, side, point):
    #     return (-1)**self._get_side(side, point) * 0.5

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

