from random import Random


class RandomMockup(Random):
    def __init__(self, value):
        super().__init__()
        self._value = value
    
    def random(self):
        return self._value
    
    @staticmethod
    def randint(a, b):
        return 0
