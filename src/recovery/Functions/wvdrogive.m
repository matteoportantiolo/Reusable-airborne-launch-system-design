function cdw = wvdrogive(diameter, length, Mach)

if Mach > 3.5       % Method not valid for higher Mach 
    Mach = 3.5;
end


sigma = 2*(180/pi)*atan(diameter/2/length);
P = (0.083 + 0.096/Mach^2)*(sigma/10)^1.69;
lod2 = (length/diameter)^2;
cdw = P*(1 -(169*lod2-16)/14/lod2/(Mach+18));

end

