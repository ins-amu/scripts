function precision=get_precision(hdr)

switch hdr.dim.datatype
    case   1,
        precision = 'ubit1';
    case   2,
        precision = 'uint8';
    case   4,
        precision = 'int16';
    case   8,
        precision = 'int32';
    case  16,
        precision = 'float32';
    case  32,
        precision = 'float32';
    case  64,
        precision = 'float64';
    case 128,
        precision = 'uint8';
    case 256 
        precision = 'int8';
    case 511 
        precision = 'float32';
    case 512 
        precision = 'uint16';
    case 768 
        precision = 'uint32';
    case 1024
        precision = 'int64';
    case 1280
        precision = 'uint64';
    case 1792,
        precision = 'float64';
    otherwise
        error('Unknown data precision.\n'); 
end

return;