select all
wireframe off
select *:P || *:T  && not hydrogen
wireframe 30
color white
select *:R && not hydrogen && sidechain
wireframe 30
color green
select *:K && not hydrogen
wireframe 50
color yellow
select *:L && not hydrogen
wireframe 50
color cpk
set picking distance
select not *:P && not *:K && not *:T && not *:L && not *:R && not hydrogen
spacefill 0.2
color cpk
select all
echo "White   = Binding Site"
echo "Green   = Rotated Side-chains"
echo "Yellow  = Known Ligand"
echo "By Atom = Docked Ligand"
