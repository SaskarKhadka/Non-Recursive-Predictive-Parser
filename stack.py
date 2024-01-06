class Stack:
    def __init__(self) -> None:
        self.__storage = []

    def values(self):
        return self.__storage

    def is_empty(self):
        return len(self.__storage) == 0

    def push(self, value):
        self.__storage.insert(0, value)

    def pop(self):
        return self.__storage.pop(0)

    def peek(self):
        return self.__storage[0]
