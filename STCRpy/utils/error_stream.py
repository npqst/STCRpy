class ErrorStream:
    def __init__(self):
        self.log = []

    def __str__(self):
        return "\n".join(self.log)

    def __repr__(self):
        return self.__str__()

    def write(self, s):
        self.log.append(str(s).strip("\n"))