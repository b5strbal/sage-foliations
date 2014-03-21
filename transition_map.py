from train_track import TrainTrack
from collections import namedtuple
from foliation import mod_one, Foliation

def new_foliation(separatrices, starting_point, starting_side,
                      is_one_sided = False, direction = 'right',
                      ending_point = None, lift_type = None):
    separatrices = sorted_separatrices(separatrices, starting_point,
                                       starting_side, is_one_sided,
                                       direction)
    if ending_point == None:
        arc_length = 1
        ending_point = starting_point
    elif direction == 'right':
        arc_length = mod_one(ending_point - starting_point)
    else:
        arc_length = mod_one(starting_point - ending_point)

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
                                          ending_point)

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
            if s1.end_side != s2.end_side:
                lengths[label] -= 1 - arc_length

    new_fol = Foliation(*gen_perm, lengths = lengths,
                     flips = flips, twist = twist)
    old_fol = separatrices[0][0].foliation
    tt_map = None if lift_type != None else get_tt_map(old_fol,
                                                       new_fol,
                                                       path_entries)
                                                       
    return (new_fol, tt_map)



def get_tt_map(old_fol, new_fol, path_entries):
    tt_new = new_fol.train_track
    tt_old = old_fol.train_track
    vertex_map = {}
    edge_map = {}

    # finding the strip which is not rectangular but L-shaped
    # more specifically, the long vertical end of it
    if not new_fol.is_bottom_side_moebius():
        long_path = (1, new_fol.num_intervals(1) - 1, 0)
    else:
        side, pos = new_fol.in_which_interval(0.5, 0)
        long_path = (0, pos, 0)
        

    for interval in new_fol.intervals():
        x1, x2 = interval.as_tuple()
        if not (x1, x2, 0) in path_entries:
            continue
        pe = [path_entries[(x1, x2, i)] for i in range(2)]
        long_end = None
        for i in range(2):
            if (x1, x2, i) == long_path:
                long_end = (i, 0)
                break
            if (pe[i].new_side, pe[i].new_end, i) == long_path:
                long_end = (i, 1)
                break
        
        a0, a1, center, b0, b1 = break_apart([pe[0].path,pe[1].path], long_end)
        v0 = vertex_map[interval] = a0[-1].end()
        v1 = vertex_map[interval.pair()] = b0[0].start()
        edge_map[tt_new.get_edge_from(interval, 'pair')] = TrainTrack.Path(center)
        edge_map[tt_new.get_edge_from(v0, 'center')] = TrainTrack.Path(a0).reversed()
        b = b1 if interval.is_flipped() else b0
        edge_map[tt_new.get_edge_from(v1, 'center')] = TrainTrack.Path(b)
        
    return TrainTrack.Map(domain = tt_new,
                          codomain = tt_old,
                          vertex_map = vertex_map,
                          edge_map = edge_map)


def break_apart(paths, long_end = None):
    diff = abs(len(paths[0]) - len(paths[1]))
    cuts = [[1, -1], [1, -1]]
    if long_end != None:
        cuts[long_end[0]][long_end[1]] += diff * (-1)**long_end[1]
    return (paths[0][:cuts[0][0]], paths[1][:cuts[1][0]],
            paths[0][cuts[0][0]:cuts[0][1]],
            paths[0][cuts[0][1]:], paths[1][cuts[1][1]:])
        
            






PathEntry = namedtuple("PathEntry", "new_side,new_i,new_end,path")

def get_pair_and_path(separatrices, side, i, end, 
                       lift_type, direction, ending_point):
    tt = separatrices[0][0].foliation.train_track
    s = separatrices[side][(i + end) % len(separatrices[side])]

    
    # on the last right end of the last interval on the top side,
    # the traintrack path has to be cut off early
    if (side, i, end) == (0, len(separatrices[0]) - 1, 1):
        endpoint = ending_point
    else:
        endpoint = None

    if direction == 'left':
        end = (end + 1) % 2
    start_end = end
    if s.is_flipped():
        end = (end + 1) % 2

    # converting first separatrix to train track path
    path = s.tt_path(end, endpoint).reversed()
    interval = path[-1].end()
    # adding connecting interval to train track path
    interval2 = interval.pair()
    path.append(tt.get_oriented_edge(interval,
                                     interval2,
                                     'pair'))
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
    s = separatrices[new_side][new_i]


    if s.is_flipped():
        end = (end + 1) % 2
    if direction == 'left':
        end = (end + 1) % 2
    if end == 1:
        new_i = (new_i - 1) % len(separatrices[new_side])

    # again, on the last right end of the last interval on the top side,
    # the traintrack path has to be cut off early
    if (new_side, new_i, end) == (0, len(separatrices[0]) - 1, 1):
        endpoint = ending_point
    else:
        endpoint = None

    # converting second separatrix to train track path
    path.extend(s.tt_path(end if direction == 'right' else (end + 1)%2, endpoint))

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
            return  mod_one(sep.endpoint - starting_point)
        else:
            return mod_one(starting_point - sep.endpoint)
        

    seps = [sorted([s for s in separatrices if 
                    s.end_side == side],
                   key = distance)
            for side in {0, 1}]

    if starting_side == 1:
        seps = reversed(seps)
    
    if is_one_sided:
        seps[0].extend(seps[1])
        seps[1] = []
        
    return seps
