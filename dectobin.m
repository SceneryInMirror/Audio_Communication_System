function bin = dectobin(dec, bit)
    bin = dec2bin(abs(dec), bit);
    for j = 1:length(dec)
        if dec(j,1) < 0
            for i = 1:bit
                if bin(j, i) == 48
                    bin(j, i) = 49;
                else
                    bin(j, i) = 48;
                end
            end
            bin(j, bit) = bin(j, bit) + 1;
            for i = bit:-1:2
                if bin(j, i) == 50
                    bin(j, i - 1) = bin(j, i - 1) + 1;
                    bin(j, i) = 48;
                end
            end
            if bin(j, 1) == 50
               bin(j, 1) = 48;
            end
        end
    end   
end