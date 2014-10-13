
class SignPerm(list):

    def __init__(self, ints):
        list.__init__(self, ints)

    def index(self, a):
        print('---')
        if a in self:
            return self.index(a)
        if -a in self:
            return self.index(-a)
        return -1


if __name__ == '__main__':

    sp = SignPerm([1, -2])
    print(type(sp), sp)
    print(sp.index(1))
    print(sp.index(-2))
    print(sp.index(2))
