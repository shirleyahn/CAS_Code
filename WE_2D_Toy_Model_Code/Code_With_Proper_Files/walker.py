class Walker:
    def __init__(self, previous_coordinates, current_coordinates, global_index, ball_center, ball_key,
                 distance_from_center, initial_step_num=0, weight=0.0, state=-1):
        self.previous_coordinates = previous_coordinates
        self.current_coordinates = current_coordinates
        self.global_index = int(global_index)  # global_index indicates the corresponding trajectory file number
        self.ball_center = ball_center
        self.ball_key = int(ball_key)
        self.distance_from_center = float(distance_from_center)
        self.initial_step_num = int(initial_step_num)  # initial_time_step indicates when the walker was born
        self.weight = float(weight)
        self.state = int(state)

    def set(self, current_coordinates, weight=0.0):
        self.current_coordinates = current_coordinates
        if weight > 0.0:
            self.weight = float(weight)

    def copy_walker(self, some_walker):
        self.previous_coordinates = some_walker.previous_coordinates
        self.current_coordinates = some_walker.current_coordinates
        self.ball_center = some_walker.ball_center
        self.ball_key = some_walker.ball_key
        self.distance_from_center = some_walker.distance_from_center
        self.initial_step_num = some_walker.initial_step_num
        self.weight = some_walker.weight
        self.state = some_walker.state
