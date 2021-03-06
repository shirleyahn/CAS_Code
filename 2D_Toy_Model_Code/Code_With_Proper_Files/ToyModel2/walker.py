class Walker:
    def __init__(self, previous_coordinates, current_coordinates, global_index, previous_ball_center,
                 current_ball_center, previous_ball_key, current_ball_key, initial_step_num=0, weight=0.0, state=-1):
        self.previous_coordinates = previous_coordinates
        self.current_coordinates = current_coordinates
        self.global_index = int(global_index)  # global_index indicates the corresponding trajectory file number
        self.previous_ball_center = previous_ball_center
        self.current_ball_center = current_ball_center
        self.previous_ball_key = int(previous_ball_key)
        self.current_ball_key = int(current_ball_key)
        self.initial_step_num = int(initial_step_num)  # initial_time_step indicates when the walker was born
        self.weight = float(weight)
        self.state = int(state)

    def __repr__(self):
        return "<Walker prev_coord:%s curr_coord:%s index:%s prev_center:%s curr_center:%s prev_key:%s " \
               "curr_key:%s step_num:%s weight:%s state:%s>" % (self.previous_coordinates,
                                                                           self.current_coordinates, self.global_index,
                                                                           self.previous_ball_center,
                                                                           self.current_ball_center,
                                                                           self.previous_ball_key,
                                                                           self.current_ball_key, self.initial_step_num,
                                                                           self.weight, self.state)

    def set(self, current_coordinates, weight=0.0):
        self.current_coordinates = current_coordinates
        if weight > 0.0:
            self.weight = float(weight)

    def copy_walker(self, some_walker):
        self.previous_coordinates = some_walker.previous_coordinates
        self.current_coordinates = some_walker.current_coordinates
        self.previous_ball_center = some_walker.previous_ball_center
        self.current_ball_center = some_walker.current_ball_center
        self.previous_ball_key = some_walker.previous_ball_key
        self.current_ball_key = some_walker.current_ball_key
        self.initial_step_num = some_walker.initial_step_num
        self.weight = some_walker.weight
        self.state = some_walker.state
