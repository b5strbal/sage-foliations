
from sage.structure.sage_object import SageObject

from matplotlib.colors import ColorConverter


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
        low_y = -0.5
        if self._foliation.is_bottom_side_moebius():
            low_y = -self.get_option('moebius_width')
        else:
            low_y = -0.5

        end_x = separatrix.intersections[-1]
        if separatrix.end_side() == 0:
            end_y = 0.5
        else:
            end_y = low_y

        draw_options = self.get_option('separatrix_draw_options')
        cc = ColorConverter()
        s = '\\definecolor{{separatrixcolor}}{{rgb}}{{{0},{1},{2}}}\n'.format(\
                *cc.to_rgb(self.get_option('separatrix_color')))
        s += '\\draw[color=separatrixcolor, {opt}] ({0},0) -- ({0},{1});\n'.format(end_x, end_y, opt = draw_options) 
        for i in range(len(separatrix.intersections) - 1):
            x = separatrix.intersections[i]
            s += '\\draw[color=separatrixcolor, {opt}] ({0},{1}) -- ({0},0.5);\n'.format(x, low_y, opt = draw_options)
        return s


#    def _tikz_of_train_track(self, train_track):
#        s = ''
#        for position in self._foliation.involution.positions():
#            s += '\\node ({name}) at ({x}, {y}) {}\n'.format(
#                    name = str(pair[0]) + '_' + str(pair[1]),
#                    x = 







    def tikz_picture(self, separatrices = []):
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

        latex.add_to_preamble('\usepackage{tikz}\n')
        s = '\\begin{{tikzpicture}}[scale = {0},'\
            'font=\\tiny]\n'.format(scale_size)

        singularities = ''
        lines = ''
        fillings = ''
        labels = ''
        if self._foliation.is_bottom_side_moebius():
            moebius_width = self.get_option('moebius_width')
            lines += '\\draw (0,-{0}) [dotted] -- (1,-{0});\n'.format(\
                    moebius_width) 
            fillings += '\\fill[yellow!{cs}!white] (0,0) rectangle '\
                    '(1,-{0});\n'.format(moebius_width, 
                            cs = color_strength / 2)
            labels += '\\node[font = \large] at (0.5, -0.1) '\
                    '{Moebius band};\n'

        for interval in self._foliation.intervals():
            if interval.side == 0:
                pos = 'above'
            else:
                pos = 'below'
            begin_percent = color_strength
            end_percent = 0

            signed_label = interval.label()
            if interval.is_flipped():
                signed_label = '-' + signed_label
                if interval > interval.pair():
                    begin_percent, end_percent = end_percent, begin_percent

            x1 = interval.endpoint(0)
            x2 = interval.endpoint(1)
            if x2 == 0:
                x2 = 1
            midx = interval.midpoint()
            y = (-1)**interval.side * 0.5
            p1 = '({0},0)'.format(x1)
            p2 = '({0},{1})'.format(x1, y)
            p3 = '({0},{1})'.format(x2, y)
            p4 = '({0},0)'.format(x2)

            color = _tikzcolor(self._foliation.index_of_label(
                interval.label()))

            lines += '\\draw {0} -- {1};\n'.format(p1, p2)
            if x1 < x2:
                lines += '\\draw[dashed] {0} -- {1};\n'.format(p2, p3)
                fillings += '\\shade[left color = {0}!{bp}!white, '\
                        'right color = {0}!{ep}!white] {1} rectangle '\
                        '{2};\n'.format(color, p1, p3, 
                                bp = begin_percent, ep = end_percent)
            else:
                lines += '\\draw[dashed] {0} -- {1};\n'.format(p2, 
                    '(1,-0.5)')
                lines += '\\draw[dashed] {0} -- {1};\n'.format(
                        '(0,-0.5)', p3)
                middle_percent = color_strength * (x2) / (1 + x2 - x1)
                fillings += '\\shade[left color = {0}!{bp}!white, '\
                        'right color = {0}!{mp}!white] {1} rectangle '\
                        '{2};\n'.format(color, p1, '(1,-0.5)', 
                                mp = middle_percent, bp = begin_percent)
                fillings += '\\shade[left color = {0}!{mp}!white, '\
                        'right color = {0}!{ep}!white] {1} rectangle '\
                        '{2};\n'.format(color, '(0,0)', p3, 
                                mp = middle_percent, ep = end_percent)

            sing_color = _tikzcolor(interval.which_singularity())
            singularities += '\\filldraw[fill={col}, draw = black] {0} '\
                    'circle (0.005);\n'.format(p2, col = sing_color)
            if length_labelling:
                labels += '\\node at ({0},0) [{1}] {{{2}}};\n'.\
                        format(midx, pos, signed_label)
            if interval_labelling:
                labels += '\\node at ({0},{1}) [{2}] {{{3}}};\n'.\
                        format(midx, y, pos, round(interval.length(), 4))

        lines += '\\draw (1,0) -- (1,0.5);\n'\
                '\\draw[very thick] (0,0) -- (1,0);\n'
        s += fillings + lines + singularities + labels
        for separatrix in separatrices:
            s += self._tikz_of_separatrix(separatrix)
        s += '\\end{tikzpicture}\n'
        return s

