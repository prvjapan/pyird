import pytest
import pyird

def test_version():
    import os
    version1 = pyird.__version__
    result = os.popen('pip show pyird').read()
    index = result.find('Version')
    version2 = result[index+9:index+14]
    #print(version1,version2)
    #print(version1==version2)
    assert version1 == version2

if __name__ == '__main__':
    test_version()
