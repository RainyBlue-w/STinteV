import re
import pkgutil

NUMBER = re.compile(r'[0-9]')
LOWER_CASE = re.compile(r'[a-z]')
UPPER_CASE = re.compile(r'[A-Z]')
OTHERS = re.compile(r'[^0-9a-zA-Z]')

# 密码强度类
class Strength:

    def __init__(self, valid, strength, message):
        self.valid = valid
        self.strength = strength
        self.message = message
    
    def __repr__(self):
        return self.strength
    
    def __str__(self):
        return self.message
    
    def __bool__(self):
        return self.valid
        
# 常用密码
def load_common_password():

    CN_common_password = pkgutil.get_data(
        __package__, 'common_passwords/passwd-CN-Top10000.txt'
    ).decode().strip().splitlines()
    EN_common_password = pkgutil.get_data(
        __package__, 'common_passwords/passwd-EN-Top10000.txt'
    ).decode().strip().splitlines()
    Keyboard_common_password = pkgutil.get_data(
        __package__, 'common_passwords/passwd-keyboard-Top500.txt'
    ).decode().strip().splitlines()
    
    return list(
        set(CN_common_password) | set(EN_common_password) | set(Keyboard_common_password)
    )

COMMON_WORDS = load_common_password()

class Password:

    TERRIBLE = 0
    SIMPLE = 1
    MEDIUM = 2
    STRONG = 3
    VERY_STRONG = 4

    @staticmethod
    def is_regular(input):
        reverse = input[::-1]
        regular = ''.join(['qwertyuiop', 'asdfghjkl', 'zxcvbnm'])
        return input in regular or reverse in regular

    @staticmethod
    def is_by_step(input):
        delta = ord(input[1]) - ord(input[0])

        for i in range(2, len(input)):
            if ord(input[i]) - ord(input[i-1]) != delta:
                return False

        return True

    @staticmethod
    def is_common(input):
        return input in COMMON_WORDS


    def __call__(self, input, min_length=8, min_types=3, level=STRONG):

        if len(input) < min_length:
            return Strength(False, 'terrible', 'Too short')

        if self.is_regular(input) or self.is_by_step(input):
            return Strength(False, 'simple', 'Regular password')

        if self.is_common(input):
            return Strength(False, 'simple', 'Common password')

        types = 0

        if NUMBER.search(input):
            types += 1

        if LOWER_CASE.search(input):
            types += 1

        if UPPER_CASE.search(input):
            types += 1

        if OTHERS.search(input):
            types += 1

        if types < 2:
            return Strength(level <= self.SIMPLE, 'simple', 'Too simple')

        if types < min_types:
            return Strength(level <= self.MEDIUM, 'medium', 'Risky password')

        if types == min_types:
            return Strength(level <= self.STRONG, 'strong', 'Strong password')
        
        return Strength(True, 'strong', 'Very strong password')

password = Password()
