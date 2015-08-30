class Walker:
    def __init__(self, coordinates, global_index, ball_center, distance_from_center, initial_step_num=0, weight=0.0,
                 state=-1):
        self.coordinates = coordinates
        self.global_index = int(global_index)  # global_index indicates the corresponding trajectory file number
        self.ball_center = ball_center
        self.distance_from_center = float(distance_from_center)
        self.initial_step_num = int(initial_step_num)  # initial_time_step indicates when the walker was born
        self.weight = float(weight)
        self.state = int(state)

    def set(self, coordinates, weight=0.0):
        self.coordinates = coordinates
        if weight > 0.0:
            self.weight = float(weight)
    
    def reset(self, ball_center, initial_step_num, weight, state):
        self.ball_center = ball_center
        self.initial_step_num = int(initial_step_num)
        self.weight = float(weight)
        self.state = int(state)

    def copy_walker(self, some_walker):
        self.coordinates = some_walker.coordinates
        self.ball_center = some_walker.ball_center
        self.distance_from_center = some_walker.distance_from_center
        self.initial_step_num = some_walker.initial_step_num
        self.weight = some_walker.weight
        self.state = some_walker.state