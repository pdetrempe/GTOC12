using AstroTime


function Epoch_to_SPICE_ET(epoch::Epoch)
    ET_J2000 = TDBEpoch(0days, origin=:j2000)
    ET = value(seconds(AstroTime.j2000(epoch))) - value(seconds(AstroTime.j2000(ET_J2000)))
end