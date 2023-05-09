# WST_kymatio
create date: 2023.05.09
by: ellen jieun oh
edited version of 'kymat.io' library to do a higher order wavelet scattering transfrom (WST)  

### original kymat.io
github link: https://github.com/kymatio/kymatio
doc link: https://www.kymat.io/

### what change?
1. scattering2d/core/scattering2d.py: can do 'order =3' ! 
2. scattering2d/filter_bank.py: only make one psi filter for each 'j' (no 'res' difference) - for simplisty 
3. scattering2d/frountend/torch_fronted.py: output shape! 
