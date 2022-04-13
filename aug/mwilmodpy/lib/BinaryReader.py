from IPython import embed
import struct

class BinaryReaderEOFException(Exception):
    def __init__(self):
        pass
    def __str__(self):
        return 'Not enough bytes in file to satisfy read request'

class BinaryReader:
    # Map well-known type names into struct format characters.
    typeNames = {
        'int8'   :'b',
        'uint8'  :'B',
        'int16'  :'h',
        'uint16' :'H',
        'int32'  :'i',
        'uint32' :'I',
        'int64'  :'q',
        'uint64' :'Q',
        'float'  :'f',
        'double' :'d',
        'char'   :'s'}

    def __init__(self, fileName):
        self.file = open(fileName, 'rb')
        
    def read(self, typeName, nElements = 1):
        typeFormat = BinaryReader.typeNames[typeName.lower()]
        typeSize = struct.calcsize(typeFormat)
        
        out= []
        for i in range(int(nElements)):
            value = self.file.read(typeSize)
            if typeSize != len(value):
                raise BinaryReaderEOFException
            out.append(struct.unpack(typeFormat, value)[0])

        if nElements == 1:
            return out[0]
        else:
            return out
    
    def close(self):
        self.file.close()

    def __del__(self):
        self.file.close()
