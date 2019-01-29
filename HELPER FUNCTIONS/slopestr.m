function out = slopestr(in)

switch in
    case '13'
        temp = '1/3 Octave';
    case '23' 
        temp = '2/3 Octave';
    case '10'
        temp = '1 Octave';
    case '00'
        temp = 'Flat Control';
    otherwise
        out = in;
        return
end % end switch

out = temp;

end % end function