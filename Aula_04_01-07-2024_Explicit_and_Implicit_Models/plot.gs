'reinit'
'set display color white'
'c'
'open ImplicitLinearAdvection1D.ctl'
*'open AdvecLinearConceitual1D.ctl'

it=0
while(it<=400)
'c'
'd fa(t='it')'
'd fc(t='it')'
if(it=399)
  'q pos'  
  it=0
endif
it=it+1
endwhile
