from __future__ import annotations
from anuga import Operator
import threading

class TestOperator(Operator):
    def __init__(self, domain, f0, fc, k, description = None, label: str | None = None, logging: bool = False):
        Operator.__init__(self, domain, description, label, logging)
        self._id: float = 0.0
        self.greet_messag = 'Hello world'
        self.__key: str = "hello"

        self._lock: threading.RLock = threading.Lock()
    
    def _generate_message(self, name: str) -> str:
        return f'Hello, {name}'
    
    def greet(self, name: str) -> str:
        return self._generate_message(name)
    
    @property
    def id(self) -> float:
        with self._lock:
            return self._id
    
    @id.setter
    def id(self, new_id: float):
        with self._lock:
            self._id = new_id

    @staticmethod
    def create() -> TestOperator:
        return TestOperator()


# test_operator

# TEST_OPERATOR

# _test_operator

# __test_operator


if __name__ == '__main__':
    a = TestOperator()
