class Const(object):
    __c0 = 299_792_458.0
    __mu0 = 1.256_637_06 * 1e-6
    __epsilon0 = 8.854_187_81 * 1e-12
    __h = 6.626_070_15 * 1e-34

    @property
    def c0(self) -> float:
        return self.__c0

    @property
    def mu0(self) -> float:
        return self.__mu0

    @property
    def epsilon0(self) -> float:
        return self.__epsilon0

    @property
    def h(self) -> float:
        return self.__h
