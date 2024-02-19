q:=578960446186580977117854925043439539266349923328202820197287920039566081674\
03;//bit_q=256
k:=2;
p:=k*q+1;

//Defining the group and the subgroup Gq 
Z:=Integers();
Zp := Integers(p);
Zp_star, iso_unit_group_to_Zp := UnitGroup(Zp);
Gq, iso_Gq_to_Zq := sub<Zp_star| ((p-1) div q)*Zp_star.1>;
length_for_tame_or_wild:=50;
tame_length:=length_for_tame_or_wild;
wild_length:=length_for_tame_or_wild;
dis_pt_length:=12;
dis_pt_length_tt:=6;
r1:=4;l1:=20;
Zq:=Integers(q);

