
from sage.structure.sage_object import SageObject

from matplotlib.colors import ColorConverter
from foliation import mod_one
from constants import TOP, BOTTOM, LEFT, RIGHT

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
            'separatrix_color':'purple',
            'separatrix_draw_options':'dashed',
            'moebius_width':0.2
            }

    def __init__(self, foliation, **options):
        self._foliation = foliation
        self._options = {}
        self._separatrices = []
        self.set_options(**options)

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

    def _tikz_of_separatrix(self, separatrix):

        end_x = separatrix.intersections[-1]
        end_y = self._end_y(separatrix.end_side())
        draw_options = self.get_option('separatrix_draw_options')
        cc = ColorConverter()
        s = '\\definecolor{{separatrixcolor}}{{rgb}}{{{0},{1},{2}}}\n'.format(\
                *cc.to_rgb(self.get_option('separatrix_color')))
        s += '\\draw[color=separatrixcolor, {opt}] ({0},0) -- ({0},{1});\n'.format(end_x, end_y, opt = draw_options) 
        for i in range(len(separatrix.intersections) - 1):
            x = separatrix.intersections[i]
            s += '\\draw[color=separatrixcolor, {opt}] ({0},{1}) --'\
                 ' ({0},0.5);\n'.format(x, self._end_y(1), opt =
                                 draw_options)
        return s

    def _end_y(self, side):
        r"""
        
        INPUT:

        - ``self`` -- 

        - ``side`` -- 

        """
        if side == 0:
            return 0.5
        if self._foliation.is_bottom_side_moebius():
            return -self.get_option('moebius_width')
        return -0.5



    def _tikz_of_train_track(self, train_track):
        s = ''
        for vertex in train_track.vertices():
            s += self._tikz_of_train_track_vertex(vertex)
        for edge in train_track.edges():
            s += self._tikz_of_train_track_edge(edge)
        return s

    def _tikz_of_train_track_vertex(self, vertex):
        interval = vertex
        return '\\fill ({x},{y}) circle (0.005);\n'.format(
            x = interval.midpoint(), y = self._end_y(interval.side)/2)



    def _tikz_of_train_track_edge(self, edge):
        s = ''
        if edge[2] == 'pair':
            s += '\\draw[->-=-0.5] ({x},{y1}) -- ({x},{y2});\n'.format(
                x = edge[0].midpoint(),
                y1 = self._end_y(edge[0].side)/2,
                y2 = self._end_y(edge[0].side))
            s += '\\draw[->-=0.5] ({x},{y1}) -- ({x},{y2});\n'.format(
                x = edge[1].midpoint(),
                y1 = self._end_y(edge[1].side),
                y2 = self._end_y(edge[1].side)/2)
            return s


        clip = edge[0].is_wrapping() or edge[1].is_wrapping()
        if clip:
            s += '\\begin{scope}\n'\
                 '\\clip (0,-0.5) rectangle (1,0.5);\n'


        m = [edge[i].midpoint() for i in {0,1}]
        y = [self._end_y(edge[i].side)/2 for i in {0,1}]
        
        shift = 0.0
        if self._foliation.is_bottom_side_moebius():
            shift = 0.5
        overlap_length = min(mod_one(edge[1].endpoint(1) + shift -
                                     edge[0].endpoint(0)),
                             mod_one(edge[0].endpoint(1) -
                                     edge[0].endpoint(0)))
        right_endpoint = mod_one(edge[0].endpoint(0) + overlap_length)

        x = []
        x.append(m[0] - edge[0].length()/2 + overlap_length/2)
        x.append(m[1] - edge[1].length()/2 +
                 mod_one(edge[0].endpoint(0) + shift -
                         edge[1].endpoint(0)) + overlap_length/2)
                             
            
        xshifts = [{0}, {0}]        
        for i in range(2):
            if x[i] < 0:
                xshifts[i].add(1)
            elif x[i] > 1:
                xshifts[i].add(-1)

        print m, x, y, overlap_length, right_endpoint
        
        s += self._center_edge_piece(m[0], y[0], x[0], 0, xshifts[0], True)
        s += self._center_edge_piece(x[1], 0, m[1], y[1], xshifts[1])

        if clip:
            s += '\\end{scope}\n'

        if self._foliation.is_bottom_side_moebius():
            color = _tikzcolor(self._foliation.index_of_label(
                edge[0].label()))

            s += '\\draw[dashed,->,{color}] ({x1},0) .. controls +(0,-0.2) '\
                 'and +(0,-0.2) .. ({x2},0);\n'.\
                 format(x1 = x[0], x2 = x[1], color = color);
        return s



    @staticmethod
    def _center_edge_piece(x1, y1, x2, y2, xshifts = {0},
                           has_arrow = False):
        r"""
        
        INPUT:

        - ``x1`` -- 

        - ``y1`` -- 

        - ``x2`` -- 

        - ``y2`` -- 

        """
        s = ''
        for xshift in xshifts:
            s += '\\draw[{arrow},xshift={shift}cm] '\
                 '({x1},{y1}) .. controls +(0,{tan}) and'\
                 ' +(0,-{tan}) .. ({x2},{y2});\n'.format(
                     x1 = x1, x2 = x2, y1 = y1, y2 = y2,
                     shift = xshift, tan = (y2 - y1)/2,
                     arrow = '->' if has_arrow else '')
        return s

        





    def tikz_picture(self, separatrices = [], train_tracks = []):
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

            x = [adjust_point(interval.endpoint(i), fol) for i in
                 [LEFT, RIGHT]]
            midx = adjust_point(interval.midpoint(), fol)
            y = [get_y(interval, i, fol) for i in [LEFT, RIGHT]]
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
            midy = (-1)**midpoint_side(interval, fol) * 0.5
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
        s += fillings + lines + singularities + labels
        for separatrix in separatrices:
            s += self._tikz_of_separatrix(separatrix)

        for train_track in train_tracks:
            s += self._tikz_of_train_track(train_track)

        s += '\\end{tikzpicture}\n'
        return s

def adjust_point(x, foliation):
    if not foliation.is_bottom_side_moebius():
        return x
    if x > 0.5:
        return 2 * (x - 0.5)
    return 2 * x

def get_side(side, point, foliation):
    if not foliation.is_bottom_side_moebius():
        return side
    if point < 0.5:
        return TOP
    else:
        return BOTTOM
    
def midpoint_side(interval, foliation):
    return get_side(interval.side, interval.midpoint(), foliation)

def get_y(interval, end, foliation):
    side = get_side(interval.side, interval.endpoint(end), foliation)
    return (-1)**side * 0.5
