from train_track import TrainTrack
from collections import namedtuple
from foliation import Foliation, SaddleConnectionError
from mymath import mod_one
from constants import *

def new_foliation(separatrices, starting_point, starting_side,
                      is_one_sided = False, direction = 'right',
                      ending_point = None, lift_type = None):
    separatrices = sorted_separatrices(separatrices, starting_point,
                                       starting_side, is_one_sided,
                                       direction)

    # print separatrices
    # print direction, starting_point, ending_point, starting_side
    if ending_point == None:
        arc_length = 1
        ending_point = starting_point
    else:
        arc_length = mod_one(ending_point - starting_point)

    flips = set()
    remaining_labels = range((len(separatrices[0]) +
                              len(separatrices[1]))/2, 0, -1)
    gen_perm = [[None] * len(separatrices[i]) for i in range(2)]
    lengths = {}
    path_entries = {}

    if gen_perm[1] == []:
        gen_perm[1] = 'moebius'
        twist = None
    else:
        twist = mod_one(separatrices[1][0].endpoint -
                        separatrices[0][0].endpoint)
        if direction == 'left':
            twist = 1 - twist

    # need only one path for each interval if it is a lift, becuase
    # then we don't need the train track map. Otherwise both
    # path_entries are needed.
    num_path_entries = 1 if lift_type != None else 2
    
    for side in range(2):
        for i in range(len(separatrices[side])):
            if gen_perm[side][i] != None:
                continue

            for end in range(num_path_entries):
                path_entries[(side,i,end)] = \
                        get_pair_and_path(separatrices,
                                          side, i, end,
                                          lift_type,
                                          direction,
                                          ending_point if direction == 'right' else
                                          starting_point,
                                          do_we_cut = arc_length < 1)

            p = path_entries[(side, i, 0)]
            label = remaining_labels.pop()
            if p.new_end == 1:
                flips.add(label)
            gen_perm[side][i] = gen_perm[p.new_side][p.new_i] = label

            s1 = separatrices[side][i]
            s2 = separatrices[side][(i+1)%len(separatrices[side])]
            if direction == 'left':
                s1, s2 = s2, s1
            lengths[label] = mod_one(s2.endpoint - s1.endpoint)
            if i == len(separatrices[side]) - 1 or s1.end_side != s2.end_side:
                # the first condition is for the case when the transverse curve
                # is orientable, the secong is for non-orientable
                lengths[label] -= 1 - arc_length

    # print lengths
    old_fol = separatrices[0][0].foliation

    x = sum(list(lengths.values()))
    y = arc_length
    if y == 1 and old_fol.is_bottom_side_moebius():
        y = 0.5
    if abs(x - y) > epsilon:
        print x, y
        print old_fol
        print lengths
        exit()

    try:
        new_fol = Foliation(*gen_perm, lengths = lengths,
                            flips = flips, twist = twist)
    except SaddleConnectionError:
        pass
        # print gen_perm, lengths, twist
        # print separatrices
        # print twist
        # # print old_fol._latex_()
        # print starting_point, ending_point
        # exit
    tt_map = None if lift_type != None else get_tt_map(old_fol,
                                                       new_fol,
                                                       path_entries)
                                                       
    return (new_fol, tt_map)



def get_tt_map(old_fol, new_fol, path_entries):
    tt_new = new_fol.train_track
    tt_old = old_fol.train_track
    vertex_map = {}
    edge_map = {}

    # Finding the strip which is not rectangular but L-shaped
    # more specifically, the long vertical end of it.
    to_be_corrected = []
    if not new_fol.is_bottom_side_moebius():
        long_path_int = (1, new_fol.num_intervals(1) - 1)
        
    else:
        side, pos = new_fol.in_which_interval(0.5, 0)
        long_path_int = (0, pos)
        
    # Also, finding the intervals on the opposite side of it that will need
    # to be corrected.
    start = new_fol.divvalues[long_path_int[0]][long_path_int[1]]
    bound = start + 0.5 if new_fol.is_bottom_side_moebius() else start
    n = new_fol.num_intervals(0) - 1
    while new_fol.divvalues[0][n] > bound:
        to_be_corrected.append((0,n))
        n -= 1


    for interval in new_fol.intervals():
        side, pos = interval.as_tuple()
        # print x1, x2
        if not (side, pos, 0) in path_entries:
            continue
        pe = [path_entries[(side, pos, i)] for i in [LEFT, RIGHT]]
        # print pe[0]
        # print pe[1]
        long_end = None
        # for i in [LEFT, RIGHT]:
        if (side, pos) == long_path_int:
            # the long end is at the beginning
            long_end = (START, LEFT)
            # print (pe[i].new_side, pe[i].new_end, i)
            # print long_path
        elif (pe[0].new_side, pe[0].new_i) == long_path_int:
            # the long end is at the end
            long_end = (END, pe[0].new_end)
        tails, center = break_apart([pe[0].path,pe[1].path], long_end)
        
        # print side, pos
        # print "tails: ", tails
        # print "center:", center
        # print "long_path_int:", long_path_int
        # print "Long end:",long_end

        # computing the path in the long part of the L-shaped strip that has to
        # be appended to some other paths on the other side
        if long_end != None:
            # finding the right tail 
            long_path_to_append = tails[long_end[0]][long_end[1]]
            
            # reversing if necessary
            if long_end[0] == END:
                long_path_to_append = TrainTrack.Path(long_path_to_append).reversed()
                
            # cutting off the first edge
            long_path_to_append = long_path_to_append[1:]

        v0 = vertex_map[interval] = tails[START][LEFT][-1].end()
        v1 = vertex_map[interval.pair()] = tails[END][LEFT][0].start()
        edge_map[tt_new.get_edge_from(interval, 'pair')] = TrainTrack.Path(center)
        edge_map[tt_new.get_edge_from(interval, 'center')] = TrainTrack.Path(tails[START][LEFT]).reversed()
        b = tails[END][RIGHT if interval.is_flipped() else LEFT]
        edge_map[tt_new.get_edge_from(interval.pair(), 'center')] = TrainTrack.Path(b)

    # print long_path_to_append
    # print to_be_corrected
    # print edge_map
    for (side, pos) in to_be_corrected:
        interval = new_fol.interval(side, pos)
        edge_map[tt_new.get_edge_from(interval, 'center')].extend(long_path_to_append)

    return TrainTrack.Map(domain = tt_new,
                          codomain = tt_old,
                          vertex_map = vertex_map,
                          edge_map = edge_map)


def break_apart(paths, long_end = None):
    # print paths[0]
    # print paths[1]
    # print long_end
    diff = abs(len(paths[LEFT]) - len(paths[RIGHT]))
    cuts = [[1, -1], [1, -1]]
    if long_end != None:
        cuts[long_end[1]][long_end[0]] += diff * (-1)**long_end[0]
    return ([[paths[0][:cuts[0][0]], paths[1][:cuts[1][0]]],
            [paths[0][cuts[0][1]:], paths[1][cuts[1][1]:]]],
            paths[0][cuts[0][0]:cuts[0][1]])

        
            






PathEntry = namedtuple("PathEntry", "new_side,new_i,new_end,path")

def get_pair_and_path(separatrices, side, i, end, 
                       lift_type, direction, ending_point, do_we_cut):
    tt = separatrices[0][0].foliation.train_track
    s0 = separatrices[side][(i + end) % len(separatrices[side])]

            
    # print separatrices
    # print side, i, end

    if direction == 'left':
        end = (end + 1) % 2
    start_end = end
    if s0.is_flipped():
        end = (end + 1) % 2

    # print s
    # print end, endpoint
    # converting first separatrix to train track path
    path = s0.tt_path(end).reversed()
    
    interval = path[-1].end()
    # adding connecting interval to train track path
    interval2 = interval.pair()
    bridge = tt.get_oriented_edge(interval, interval2, 'pair')
    path.append(bridge)
    interval = interval2
    if interval.is_flipped():
        end = (end + 1) % 2
    if end == 1:
        interval = interval.next()
    is_flipped_so_far = start_end != end 
    new_side, new_i = matching_sep_index(separatrices,
                                              interval,
                                              lift_type,
                                              is_flipped_so_far, side)
    s1 = separatrices[new_side][new_i]
    
    if s1.is_flipped():
        end = (end + 1) % 2
    if direction == 'left':
        end = (end + 1) % 2
    if end == 1:
        new_i = (new_i - 1) % len(separatrices[new_side])

    # print new_side, new_i, end
    # print ending_point


    # converting second separatrix to train track path
    path.extend(s1.tt_path(end if direction == 'right' else (end + 1)%2))


    # we need to collect intersections for the special case when we need to
    # cut something off from a path

    # on the last right end of the last interval on the top side,
    # the traintrack path has to be cut off early
    begin_cut = (side, i, end) == (0, len(separatrices[0]) - 1, 1)
    end_cut = (new_side, new_i, end) == (0, len(separatrices[0]) - 1, 1)

    # do_we_cut: when the transformation is a rotation or reversing, we
    # should not cut
    if not (begin_cut or end_cut) or not do_we_cut:
        return PathEntry(new_side,new_i,end,path)

    # print '------------'
    # print path
    # print '------------'

    p0 = p1 = None
    if ending_point in s0.intersections():
        if end_cut:
            path = s0.tt_path(end, ending_point, 'after').reversed()
            return PathEntry(new_side,new_i,end,path)
        p0 = s0.tt_path(end, ending_point, 'before').reversed()
    else:
        if begin_cut:
            path = s1.tt_path(end, ending_point, 'after')
            return PathEntry(new_side,new_i,end,path)
        p1 = s1.tt_path(end, ending_point, 'before')

    if p0 == None:
        p0 = s0.tt_path(end).reversed()
    if p1 == None:
        p1 = s1.tt_path(end)
        
    path = p0 + [bridge] + p1

    return PathEntry(new_side,new_i,end,path)



def matching_sep_index(separatrices, interval, lift_type,
                       is_flipped_so_far, orig_side):
    for side in range(2):
        for j in range(len(separatrices[side])):
            s = separatrices[side][j]
            is_total_flipped = (s.is_flipped() != is_flipped_so_far)
            if s.first_interval(0) == interval:
                
                if lift_type == None or \
                   lift_type == 'foliation' and not is_total_flipped or \
                   lift_type == 'surface' and \
                   (s.end_side == orig_side) == is_total_flipped:
                    return (side, j)
    
    assert(False)



def sorted_separatrices(separatrices, starting_point, starting_side,
                        is_one_sided, direction = 'right'):
    def distance(sep):
        if direction == 'right':
            return mod_one(sep.endpoint - starting_point)
        else:
            return mod_one(starting_point - sep.endpoint)
        

    seps = [sorted([s for s in separatrices if 
                    s.end_side == side],
                   key = distance)
            for side in {0, 1}]

    if starting_side == 1:
        seps = list(reversed(seps))
    
    if is_one_sided:
        seps[0].extend(seps[1])
        seps[1] = []
        
    return seps
