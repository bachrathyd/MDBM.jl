module MDBM


export mdbm_problem, mdbm_object, refine!, checkncube!, checkneighbour!, interpolate!, DTconnect!

function axdoubling!(ax::Array{Array{Float64,1}},dims=[0])#::Int64)
dimindex=dims[1]==0 ? collect(1:length(ax)) : dims
  for i in dimindex
    ax[i]=sort(vcat(ax[i],ax[i][1:end-1]+diff(ax[i])/2));
    #println(length(ax[i]))
  end
end

function axdoubling!(ax,dims=[0])
  dimindex=dims[1]==0 ? collect(1:length(ax)) : dims
    for i in dimindex
    ax[i]=sort(vcat(ax[i],ax[i][1:end-1]+diff(ax[i])/2));
    #println(length(ax[i]))
  end

end

function indtosub(Nax::Array{Int64,1},linindex::Array{Int64,1})
  if isempty(linindex)
     return zeros(Int64,length(Nax),0)
  else
    Ndim=length(Nax)
    NaxNprod=[1;cumprod(Nax[1:end-1])]
    A=eye(Int64,Ndim)+diagm(-Nax[1:(Ndim-1)],1)
    fullindex=ceil.(Int64,((1 ./NaxNprod)*transpose(linindex))-(1/prod(Nax))*0.25)
    #1e-12 is for fremove the numerical error in case of ineger numbers!!!
    #indexvect=(A*(fullindex-1)+1)
    return (A*(fullindex-1)+1)
  end
end

function subtoind(Nax::Array{Int64,1},vectindex::Array{Int64,2})
  if isempty(vectindex)
    return zeros(Int64,0)
  else
    NaxNprod=[1;cumprod(Nax[1:end-1])]
    return ((transpose(vectindex)-1)*NaxNprod)+1
  end
end

function subtoind(Nax::Array{Int64,1},vectindex::Array{Int64,3})
  if isempty(vectindex)
    return zeros(Int64,0)
  else
    #  NaxNprod=[1;cumprod(Nax[1:end-1])].'
    #return  mapslices(x->(NaxNprod*(x-1))+1,vectindex,[1,2])
    NaxNprod=[1;cumprod(Nax[1:end-1])]

#=
"@devec"
    @time begin
     resout = zeros(1,size(vectindex,2),size(vectindex,3))
     for m = 1 : size(vectindex,3)
     for j = 1 : size(vectindex,2)
       for i = 1 :  length(NaxNprod)
         resout[1,j,m] += (vectindex[i,j,m] -1)* NaxNprod[i] + 1
       end
     end
   end
   end
=#

    return  sum((vectindex-1).*NaxNprod,1)+1
  end
end
function subtoind!(res::Array{Int64,2},Nax::Array{Int64,1},vectindex::Array{Int64,2})
  #res=Array{Int64}(size(vectindex,2))
  if isempty(vectindex)
    return zeros(Int64,0)
  else
    NaxNprod=[1;cumprod(Nax[1:end-1])]
    res=((transpose(vectindex)-1)*NaxNprod)+1 # TODO: melyik a gyorsabb?
    #res=sum((vectindex-1).*NaxNprod,1)+1 # TODO: melyik a gyorsabb?
  end
end

function subtoind(Nax::Array{Int64,1},vectindex::Array{Int64,1})#copmputing for one element only (vectindex is a vector and not a matrix)
  if isempty(vectindex)
    return zeros(Int64,0)
  else
    NaxNprod=[1;cumprod(Nax[1:end-1])]
    out=Int64
    return out=(((transpose(vectindex)-1)*NaxNprod)+1)
  end
end

function vect2pos(ax::Array{Array{Float64,1}},vectindex::Array{Int64,2})
  [ax[i][vectindex[i,j]] for i= 1:length(ax) , j=1:size(vectindex,2)]
end
function vect2pos(ax::Array{Array{Float64,1}},vectindex::Array{Int64,1})
  [ax[i][vectindex[i]] for i= 1:length(ax)]
end
function vect2pos(ax::Array{Array{Float64,1}},vectindex::Array{Int64,2},P::Array{Float64,2})#linearly interpolates within the cube
  #0<P<1
  [ax[i][vectindex[i,j]]+P[i,j]*(ax[i][vectindex[i,j+1]]-ax[i][vectindex[i,j]]) for i= 1:length(ax) , j=1:size(vectindex,2)]
end
function vect2pos(ax::Array{Array{Float64,1}},vectindex::Array{Int64,1},P::Array{Float64,1})#linearly interpolates within the cube
  #0<P<1
  [ax[i][vectindex[i]]+P[i]*(ax[i][vectindex[i]+1]-ax[i][vectindex[i]]) for i= 1:length(ax)]
end

function vect2pos(ax::Array{Array{Float64,1}},vectposition::Array{Float64,2})#linearly interpolates within the cube
  #(vectposition[i,j]%1) #relative index within the cube
  #vectindex=floor(Int64,vectposition[i,j])
  [ax[i][floor(Int64,vectposition[i,j])]*(1-(vectposition[i,j]%1))+ax[i][floor(Int64,vectposition[i,j])+1]*(vectposition[i,j]%1) for i= 1:length(ax) , j=1:size(vectindex,2)]
end
function vect2pos(ax::Array{Array{Float64,1}},vectposition::Array{Float64,1})#linearly interpolates within the cube
  #(vectposition[i,j]%1) #relative index within the cube
  #vectindex=floor(Int64,vectposition[i,j])
  [ax[i][floor(Int64,vectposition[i])]*(1-(vectposition[i]%1))+ax[i][floor(Int64,vectposition[i])+1]*(vectposition[i]%1) for i= 1:length(ax)]
end


function db_combnk(v,k)
m = size(v[:],1);
if m==0
    c=[];
    #     disp('the vector is empty')
else
    if k<1
        #         disp('positive integer element has to extract')
        c=[];
            println("ez lefut?");
      #      prinlt(c)
    elseif k==1
        c=v[:];
    else
        if k>m
            c=[];
          #            disp('vector has less elements than needed')
        else
            cnew=Int64[];
            for kallc=1:m
                cadditional=db_combnk(collect(kallc+1:m),k-1);
                if ~isempty(cadditional)
                  cnewpart=Int64[kallc*ones(size(cadditional,1),1) cadditional];
                  if isempty(cnew)
                    cnew=cnewpart
                  else
                    cnew=vcat(cnew,cnewpart);
                  end
                end
            end
            c=cnew;

            c=v[c];
            if k==1
                c=c';
            end
        end
    end
end
return c;
end


type mdbm_object
  f::Function## evaluated function
  fconstrain::Function## evaluated function
  fvectorized::Bool ## is function f can be called in a vectorized form
  ax::Array{Array{Float64,1},1}##description of the initial grid
  Ndim::Int64
  Nax::Array{Int64,1}
  Naxstride::Array{Int64,1}

  Ncodim::Int64
  HC::Array{Float64,2}#T computed function values and constrain value
  linindex::Array{Int64,1} # corresponding linear-indices
  vectindex::Array{Int64,2}# corresponding sub-indices
  #N::Int64
  #compind
  #pointerp
  #DT
  ncubelin::Array{Int64,1} # linear-indices of the bracketing n-cubes
  ncubevect::Array{Int64,2} # sub-indices of the bracketing n-cubes

  posinterp::Array{Float64,2}# the interpolated valuse within the bracketing n-cubes
  gradient::Array{Float64,3}# the corresponding gradients within the bracketing n-cubes

  DT1::Array{Int64,2}#line 'tiangulation' of the resultant interpolated values (basd on the n-cube sub-indices)
  DT2::Array{Int64,2}#surface tiangulation of the resultant interpolated values (basd on the n-cube sub-indices)

  isconstrained::Bool #is f provides contrain? all the constraints are combinded!!! length C===1
  interporder::Int64 # 0,1,2
  selectionmode::Int64 # 0-safe selection, 1st order interpolation mased,2???
end


function mdbm_problem(f::Function,axabstract:: Array{AbstractArray{T,1} where T,1};fconstrain::Function=(x...)->1.0,isconstrained::Bool=true,Ncodim::Int64=0,interporder::Int64=1,selectionmode::Int64=0, fullrefinenum::Int64=0, fvectorized=false)
  ax=[ convert(Array{Float64,1},collect(axloc))  for axloc in axabstract]

  Ndim=length(ax)
  Nax=[length(kax) for kax in ax]
  Naxstride=[1;cumprod(Nax[1:end-1])]

  linindex=collect(1:prod([length(kax) for kax in ax]))
  vectindex=indtosub(Nax,linindex)
  if Ncodim==0
    if fvectorized
      ax1loc=[axloc[1] for axloc in ax];
      ax1loc=reshape(ax1loc, Ndim, 1)
      Ncodim=length([f(ax1loc...);fconstrain(ax1loc...)])-isconstrained#TODO: paraméterátadás, ha kell több bemenet
    else
      Ncodim=length([f([axloc[1] for axloc in ax]...);fconstrain([axloc[1] for axloc in ax]...)])-isconstrained#TODO: paraméterátadás, ha kell több bemenet
    end
  end

  #HC=Array{Float64}(Ncodim+isconstrained,prod([length(kax) for kax in ax]))
  #HC=mapslices(f,vect2pos(ax,vectindex),[1])
  HC=fevaluate(f,fconstrain,ax,vectindex,Ncodim+isconstrained,isfvectorized=fvectorized)

  ncubelin=collect(1:prod([length(kax)-1 for kax in ax]))
  ncubevect=indtosub(Nax-1,ncubelin)
  #N::Int64
  #compind
  #pointerp
  #DT
  posinterp=Array{Float64}(length(Nax),0)
  gradient=Array{Float64}(Ncodim+isconstrained,length(Nax),0)
  DT1=Array{Int64}(2,0)
  DT2=Array{Int64}(3,0)
  mdbm=mdbm_object(f,fconstrain,fvectorized,deepcopy(ax),Ndim,Nax,Naxstride,Ncodim,HC,linindex,vectindex,ncubelin,ncubevect,posinterp,gradient,DT1,DT2,isconstrained,interporder,selectionmode)

  if fullrefinenum>0
    checkncube!(mdbm)
    #interpolate!(mdbm)
    for k=1:fullrefinenum
      refine!(mdbm)
      checkncube!(mdbm)
    end
    checkneighbour!(mdbm)
    interpolate!(mdbm)
  end
  return mdbm
end


function fevaluate(f::Function,fconstrain::Function,ax::Array{Array{Float64,1},1},vectindex::Array{Int64,2},NcodimNplusisconstrained=1;isfvectorized=false)
  #compute the elements even if it is already computed before
  #check and filter the vectindex to elimiate the redundant computation
  if size(vectindex,2)==0
    HC=zeros(Float64,NcodimNplusisconstrained,0)
  else
    axialpos=vect2pos(ax,vectindex)
    HC=zeros(Float64,NcodimNplusisconstrained,size(vectindex,2))
    if isfvectorized
      HC=[f(axialpos);fconstrain(axialpos)]
    else
      @simd for k=1:size(vectindex,2)
      HC[:,k]=[f(axialpos[:,k]...);fconstrain(axialpos[:,k]...)]
      end
    end
   #  @time HC=mapslices(f,vect2pos(ax,vectindex),[1]) # slower
  end
  return HC
end

function isbracketing(Nax::Array{Int64,1},HC::Array{Float64,2},isconstrained::Bool=false,selectionmode::Int64=0)
  #TODO: create the higher order selection mode
  #HCsign=HC[1:end-isconstrained,:].>0
  if isconstrained
    HsignPOS=HC[1:end-isconstrained,:].>=0
    HsignNEG=HC[1:end-isconstrained,:].>=0
    CsignPOS=HC[end,:].>=0
    return (all(any(HsignPOS,2).&any(HsignNEG,2))).&(any(CsignPOS[:]))
  else
    HsignPOS=HC[1:end,:].>=0
    HsignNEG=HC[1:end,:].<=0
    return all(any(HsignPOS,2).&any(HsignNEG,2))
  end
end

function isbracketing(Nax::Array{Int64,1},HC::Array{Float64,3},isconstrained::Bool=false,selectionmode::Int64=0)
  #TODO: create the higher order selection mode
  #HCsign=HC[1:end-isconstrained,:].>0
  if isconstrained
    HsignPOS=HC[1:end-isconstrained,:,:].>=0
    HsignNEG=HC[1:end-isconstrained,:,:].<=0
    CsignPOS=HC[end,:,:].>=0#dimensions reducred, becuse it is [1xnxm]->[nxm]  ## if all is negative, then it is a false ncube
    #@time (all(any(Hsign,2).&any(.!Hsign,2),1))[:].&(any(Csign,1))[:]
    #@time squeeze(squeeze(all(any(Hsign,2).&any(.!Hsign,2),1),1),1).&squeeze(any(Csign,1),1)
    return squeeze(squeeze(all(any(HsignPOS,2).&any(HsignNEG,2),1),1),1).&squeeze(any(CsignPOS,1),1)
  else
    HsignPOS=HC[1:end,:,:].>=0
    HsignNEG=HC[1:end,:,:].>=0
    #@time all(any(Hsign,2).&any(.!Hsign,2),1)[:]
    #@time squeeze(squeeze(all(any(Hsign,2).&any(.!Hsign,2),1),1),1)
    return squeeze(squeeze(all(any(HsignPOS,2).&any(HsignNEG,2),1),1),1)
  end
end

#=
function isbracketing(mdbm::mdbm_object)
  #TODO: create the higher order selection mode
  return isbracketing(mdbm.Nax,mdbm.HC,mdbm.isconstrained,mdbm.selectionmode)
end
=#
function checkncube!(mdbm::mdbm_object)
  propercube=checkncube(mdbm::mdbm_object,mdbm.ncubevect)
  mdbm.ncubelin=mdbm.ncubelin[propercube]
  mdbm.ncubevect=mdbm.ncubevect[:,propercube]
end

function checkncube(mdbm::mdbm_object,ncubevect::Array{Int64,2})
  if isempty(ncubevect)
    println("There is no bracketing n-cubes to check!")
    propercube=Array{Bool}(0)
  else
    sub2=indtosub(2*ones(Int64,mdbm.Ndim),collect(1:2^mdbm.Ndim))-1
    sub22=Array{Int64}(size(sub2,1),1,size(sub2,2))
    sub22[:]=sub2[:]
    cubinds=permutedims(subtoind(mdbm.Nax,ncubevect.+sub22),[3,2,1])

    locind=Array{Int64}(2^mdbm.Ndim)
    propercube=Array{Bool}(size(ncubevect,2))

    #@time cubindsU=unique(cubinds)
    p = sortperm(cubinds[:])
    pp = sortperm(p)
    cubindsU=unique(cubinds[p])
    #cubinds[p]
    unixind=indexin_sorted(cubinds[p],cubindsU)
    #cubindsU[unixind[pp]]

    HCindU=indexin_sorted(cubindsU,mdbm.linindex)

    HCallcubes=mdbm.HC[:,HCindU[unixind[pp]]]
    HCallcubes=reshape(HCallcubes,mdbm.Ncodim+mdbm.isconstrained,2^mdbm.Ndim,size(ncubevect,2))
    #@time propercube=[isbracketing(mdbm.Nax,HCallcubes[:,:,j],mdbm.isconstrained,mdbm.selectionmode) for j=1:size(ncubevect,2)]
    #@time propercube=mapslices(X ->isbracketing(mdbm.Nax,X,mdbm.isconstrained,mdbm.selectionmode), HCallcubes, [1,2])

    propercube=isbracketing(mdbm.Nax,HCallcubes,mdbm.isconstrained,mdbm.selectionmode)


    #cubeindloc=Array{Int64}(2^mdbm.Ndim)
    #HC=Array{Int64}(mdbm.Ncodim+mdbm.isconstrained,2^mdbm.Ndim)

    #@simd for j=1:size(ncubevect,2)
    #  propercube[j]=isbracketing(mdbm.Nax,mdbm.HC[:,indexin_sorted(cubinds[:,j,1],mdbm.linindex)],mdbm.isconstrained,mdbm.selectionmode)
    #end
    if isempty(mdbm.ncubelin)
      println("There is no bracketing n-cubes!")
    end
  end
  return propercube
end

function refine!(mdbm::mdbm_object,dims=[0])
   dimindex=dims[1]==0 ? collect(1:length(mdbm.ax)) : dims
   axdoubling!(mdbm.ax,dimindex)
   mdbm.Nax[dimindex]=(mdbm.Nax[dimindex]-1)*2+1
   mdbm.Naxstride=[1;cumprod(mdbm.Nax[1:end-1])]

   mdbm.vectindex[dimindex,:]=(mdbm.vectindex[dimindex,:]-1)*2+1
   mdbm.linindex=subtoind(mdbm.Nax,mdbm.vectindex)

   mdbm.ncubevect[dimindex,:]=(mdbm.ncubevect[dimindex,:]-1)*2+1
  # expanding along each selected dimension
   @simd for kdim in dimindex
    mdbm.ncubevect=hcat(mdbm.ncubevect,mdbm.ncubevect.+(collect(1:mdbm.Ndim).==kdim))
  end
  mdbm.ncubelin=subtoind(mdbm.Nax,mdbm.ncubevect)
  #sorting

  p = sortperm(mdbm.ncubelin)
  mdbm.ncubelin=mdbm.ncubelin[p]
  mdbm.ncubevect=mdbm.ncubevect[:,p]

  evalmissingnodes!(mdbm)

end

function evalmissingnodes!(mdbm::mdbm_object)
  evalmissingnodes!(mdbm::mdbm_object,mdbm.ncubevect)
end

function evalmissingnodes!(mdbm::mdbm_object,ncubevect::Array{Int64,2})
  if isempty(ncubevect)
    println("There is no bracketing n-cubes!")
  else
    sub2=indtosub(2*ones(Int64,mdbm.Ndim),collect(1:2^mdbm.Ndim))-1
    sub22=Array{Int64}(size(sub2,1),1,size(sub2,2))
    sub22[:]=sub2[:]
    cubinds=permutedims(subtoind(mdbm.Nax,ncubevect.+sub22),[3,2,1])

    cubinds=unique(cubinds[:])
    sort!(cubinds)
    #bbbbbbbbbbbbbbb
    #newlocind=cubinds[indexin(cubinds,mdbm.linindex).==0]

    newlocind=cubinds[indexin_sorted_iszero(cubinds,mdbm.linindex)]
    newvectind=indtosub(mdbm.Nax,newlocind)
    mdbm.linindex=vcat(mdbm.linindex,newlocind)
    mdbm.vectindex=hcat(mdbm.vectindex,newvectind)
    if size(newvectind,2)==1
      newvectind=reshape(newvectind, mdbm.Ndim, 1) # force it to a matrix type
    end

    mdbm.HC=hcat(mdbm.HC,fevaluate(mdbm.f,mdbm.fconstrain,mdbm.ax,newvectind,mdbm.Ncodim+mdbm.isconstrained,isfvectorized=mdbm.fvectorized))
    #TODO: megnézni a sorbarakást, hogy kell-e, gyorsabb lesz-e
    #sorting
    p = sortperm(mdbm.linindex)
    mdbm.HC=mdbm.HC[:,p]
    mdbm.linindex=mdbm.linindex[p]
    mdbm.vectindex=mdbm.vectindex[:,p]
  end
end

#TODO: @inbounds ????
function interpolate!(mdbm::mdbm_object,bracketingoverwrite::Bool=true)
if isempty(mdbm.ncubelin)
  println("There is no bracketing n-cubes!")
  return
end
propercube=trues(size(mdbm.ncubevect,2))
  if mdbm.interporder==0
  #println("zeroth order interp")
  mdbm.posinterp=Array{Float64}(mdbm.Ndim,size(mdbm.ncubevect,2))
  mdbm.gradient=Array{Float64}(mdbm.Ncodim,mdbm.Ndim,0)
  sub2=indtosub(2*ones(Int64,mdbm.Ndim),collect(1:2^mdbm.Ndim))-1
  @simd for j=1:size(mdbm.ncubevect,2)
    #mdbm.posinterp[:,j]=vect2pos(ax,mdbm.ncubevect[:,j]) #"lower left corener"
    #mdbm.posinterp[:,j]=median(vect2pos(ax,mdbm.ncubevect[:,j].+sub2),2) #midpoint
    mdbm.posinterp[:,j]=vect2pos(mdbm.ax,mdbm.ncubevect[:,j],0.5*ones(mdbm.Ndim))  #much faster!!!!
  end

elseif mdbm.interporder==1
  #println("first order interp")
  mdbm.posinterp=Array{Float64}(mdbm.Ndim,size(mdbm.ncubevect,2))
  mdbm.gradient=Array{Float64}(mdbm.Ndim,mdbm.Ncodim,size(mdbm.ncubevect,2))

  sub2=indtosub(2*ones(Int64,mdbm.Ndim),collect(1:2^mdbm.Ndim))-1
  sub22=Array{Int64}(size(sub2,1),1,size(sub2,2))
  sub22[:]=sub2[:]
  cubinds=permutedims(subtoind(mdbm.Nax,mdbm.ncubevect.+sub22),[3,2,1])
  #sub2=indtosub(2*ones(Int64,mdbm.Ndim),collect(1:2^mdbm.Ndim))-1
  #cubinds=Array{Int64}(2^mdbm.Ndim)
  locind=Array{Int64}(2^mdbm.Ndim)

  axlocdimles=Array{Float64}(0,1)
    for kdim=1:mdbm.Ndim
        axlocdimles=hcat([axlocdimles;-ones(1,2^(kdim-1))],[axlocdimles;ones(1,2^(kdim-1))]);
    end
    TAn=hcat(-ones(2^mdbm.Ndim,1),transpose(axlocdimles))
    #     TAn2=inv(TAn.'*TAn);
    #     TAtrafo=TAn2*TAn.';
    TAtrafo=(transpose(TAn)*TAn)\transpose(TAn)


As=Array{Float64}(mdbm.Ncodim)
ns=Array{Float64}(mdbm.Ndim,mdbm.Ncodim)
solall=Array{Float64}(mdbm.Ncodim+mdbm.Ndim,mdbm.Ncodim)
ConstrainDominant=Bool
HCloc=Array{Float64}(mdbm.Ndim+mdbm.isconstrained,2^mdbm.Ndim)
P=Array{Float64}(mdbm.Ndim)
  #println("1sr oredr interpolation time!")
#@simd

for j=1:size(mdbm.ncubevect,2)
  locind=indexin_sorted(cubinds[:,j,1],mdbm.linindex)
  #  HC=mdbm.HC[:,locind]
  #the functionvalue of the constraint must be considered during the interpolation to force the point into the boundary of the surface
  HCloc=mdbm.HC[:,locind]
  ConstrainDominant=(mdbm.isconstrained && any(HCloc[end,:].<=0));
  solall=TAtrafo*transpose(HCloc[1:mdbm.Ncodim+ConstrainDominant,:]);

  As=solall[1,:];#it is not a real distance within the n-cube (it is ~n*A)!!!
  ns=solall[2:end,:]
  P=ns*((transpose(ns)*ns)\As);
  #TODO: what if it falls outside of the n-cube
  #TODO: it should be removed ->what shall I do with the bracketing cubes?
            if norm(P,20)>2 #  bondcube(kcubes)=all(norm(P)<=(Ndim^0.5));
                  propercube[j]=false;# P=2*P/norm(P,20);
            end
  mdbm.gradient[:,:,j]=ns[1:mdbm.Ndim,1:mdbm.Ncodim]
  #mdbm.posinterp[:,j]=vect2pos(mdbm.ax,mdbm.ncubevect[:,j]+(P/2+0.5)) Can provide points outside of the range!!!
  #mdbm.posinterp[:,j]=vect2pos(mdbm.ax,min(mdbm.Nax-1e-8,max(ones(mdbm.Ndim),mdbm.ncubevect[:,j]+(P/2+0.5))))#force it within the cube
  mdbm.posinterp[:,j]=vect2pos(mdbm.ax,mdbm.ncubevect[:,j],(P/2+0.5)) # if P falls outside the cube, only the gradient is used to extrapolate!!!
  end
else #mdbm.interporder==2
end
if bracketingoverwrite
mdbm.gradient=mdbm.gradient[:,:,propercube]
mdbm.posinterp=mdbm.posinterp[:,propercube]
mdbm.ncubevect=mdbm.ncubevect[:,propercube]
mdbm.ncubelin=mdbm.ncubelin[propercube]
end
end

function indexin_sorted(a::Array{Int64,1},b::Array{Int64,1})
  # b must contain all the lements of a
  # a and b must be sorted
  if isempty(a)
    return zeros(0)
  elseif length(b)==1 ##much faster!
    return (a.==b)*1
  else
  out=zeros(Int64,size(a))
  leng::Int64=length(b)

  q::Int64=1
  q1::Int64=1
  q2::Int64=length(b)
  for k= 1:length(a)

    if (q2-q1)==1
      q=(b[q1]==a[k]) ? q1 : q2
    else
    while (b[q]!=a[k]) & (q2>q1+1)
      if b[q]>a[k]
        q2=q
        q=floor(Int64,(q+q1)/2)
      else
        q1=q
        q=ceil(Int64,(q+q2)/2)
      end
      #print([q1;q;q2])
      #print([b[q]])
      #println([a[k]])
    end
  end
    out[k]=b[q]==a[k]?q:0
    q=max(out[k],q1)
    q1=q
    q2=length(b)
    #print("-----")
    #print(out[1:k])
    #println("----")
    if q1>length(b)#all the element is larger than the last one
      break
    end
  end

  return out
  end
end

function indexin_sorted_iszero(a::Array{Int64,1},b::Array{Int64,1})
  # b must contain all the lements of a
  # a and b must be sorted
  if isempty(a)
    return zeros(0)
  elseif length(b)==1 ##much faster!
    return (a.==b)*1
  else
  out=falses(size(a))
  leng::Int64=length(b)

  q::Int64=1
  q1::Int64=1
  q2::Int64=length(b)
  for k= 1:length(a)

    if (q2-q1)==1
      q=(b[q1]==a[k]) ? q1 : q2
    else
    while (b[q]!=a[k]) & (q2>q1+1)
      if b[q]>a[k]
        q2=q
        q=floor(Int64,(q+q1)/2)
      else
        q1=q
        q=ceil(Int64,(q+q2)/2)
      end
      #print([q1;q;q2])
      #print([b[q]])
      #println([a[k]])
    end
  end
    out[k]=b[q]==a[k] ? false : true
    q=max(out[k],q1)
    q1=q
    q2=length(b)
    #print("-----")
    #print(out[1:k])
    #println("----")
    if q1>length(b)#all the element is larger than the last one
      break
    end
  end

  return out
  end
end

function checkneighbour!(mdbm::mdbm_object)
  if isempty(mdbm.ncubelin)
    println("There is no bracketing n-cubes to check!")
  else
    sub3=indtosub(3*ones(Int64,mdbm.Ndim),collect(1:3^mdbm.Ndim))-2
    sub33=Array{Int64}(size(sub3,1),1,size(sub3,2))
    sub33[:]=sub3[:]
    newbracketinncubes=mdbm.ncubelin
    checkedncubes=mdbm.ncubelin
    bracketingcube=true
    while any(bracketingcube)
      subvals=indtosub(mdbm.Nax,newbracketinncubes).+sub33
      #subvals=subvals[:,:]
      subvals=reshape(subvals,size(subvals,1),:)
      porperrange=all(vcat(subvals.>(mdbm.Nax*0),subvals.<(mdbm.Nax)),1) #TODO: >0 vagy >=0
      #subvals=subvals[:,porperrange[:]]
      cubinds=subtoind(mdbm.Nax,subvals[:,porperrange[:]])
      sort!(cubinds)
      cubinds=unique(cubinds[:]) #these neighbours must be checked
      #issorted(cubinds)
      #issorted(checkedncubes)
      newneighbour=cubinds[indexin_sorted_iszero(cubinds,checkedncubes)]
      if !isempty(newneighbour)
        evalmissingnodes!(mdbm,indtosub(mdbm.Nax,newneighbour))
        #issorted(newneighbour)
        bracketingcube=checkncube(mdbm,indtosub(mdbm.Nax,newneighbour))
        newbracketinncubes=newneighbour[bracketingcube]
        checkedncubes=sort(vcat(newneighbour,checkedncubes))
        mdbm.ncubelin=vcat(mdbm.ncubelin,newbracketinncubes)
      else
        bracketingcube=false
      end
    end
    sort!(mdbm.ncubelin)
    mdbm.ncubevect=indtosub(mdbm.Nax,mdbm.ncubelin)
  end
end

function DTconnect!(mdbm::mdbm_object)
  ncubevect=mdbm.ncubevect
  ncubelin=mdbm.ncubelin
if (mdbm.Ndim-mdbm.Ncodim)>0
  #TODO: compute the higher order conncetion (for triangles "szomszédossági mátrix gárffal jó lesz!!!)")

  closeneigh2=diagm(ones(Int64,mdbm.Ndim))[:,:]
  #closeneigh2=cat(3,diagm(ones(Int64,mdbm.Ndim)),-diagm(ones(Int64,mdbm.Ndim)))[:,:]#TODO: "left neighbour is not neseccasry, this connection will be detected when the "left neighbout" is analised
  closeneigh22=Array{Int64}(size(closeneigh2,1),1,size(closeneigh2,2))
  closeneigh22[:]=closeneigh2[:]
  cubinds=permutedims(subtoind(mdbm.Nax,ncubevect.+closeneigh22),[3,2,1])
  DT1=Array{Int64,2}(0,2)
   for j=1:size(ncubevect,2) #@simd
    locind=indexin_sorted(cubinds[:,j,1],ncubelin)
    #DT1=vcat(DT1,hcat(ncubelin[j]*ones(Int64,sum(locind.!=0)),ncubelin[locind[locind.!=0]]))
    DT1=vcat(DT1,hcat(j*ones(Int64,sum(locind.!=0)),locind[locind.!=0]))
  end
  mdbm.DT1=DT1
end

if (mdbm.Ndim-mdbm.Ncodim)>1 #at least surface

  cornerneigh2=(indtosub(2*ones(Int64,mdbm.Ndim),collect(1:2^mdbm.Ndim))-1)[:,2:end]#TODO: -1:1 is necessary, it is possible to lose some connections!!!, but if -1:1 is used, then every squares will be plotted by 4 triangles!!!
  cornerneigh22=Array{Int64}(size(cornerneigh2,1),1,size(cornerneigh2,2))
  cornerneigh22[:]=cornerneigh2[:]
  cubinds=permutedims(subtoind(mdbm.Nax,ncubevect.+cornerneigh22),[3,2,1])
  DT1corner=Array{Int64,2}(0,2)
  @simd for j=1:size(ncubevect,2)
    locind=indexin_sorted(cubinds[:,j,1],ncubelin)
    #DT1=vcat(DT1,hcat(ncubelin[j]*ones(Int64,sum(locind.!=0)),ncubelin[locind[locind.!=0]]))
    DT1corner=vcat(DT1corner,hcat(j*ones(Int64,sum(locind.!=0)),locind[locind.!=0]))
  end
  DT1corner
  #mdbm.DT1=DT1corner

  #i=[1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8]
  #j=[2,3,4,1,3,6,1,2,4,1,3,5,4,6,7,2,5,7,5,6,8,7]
  #v=ones(Int64,size(j))
  #A=full(sparse(i,j,v))
  #A=(sparse(i,j,v))

  #DT-elő kell állítani, sarok menti szomszédokkal is és akkor már jó is lesz!!!


  A=sparse(DT1corner[:,1],DT1corner[:,2], ones(Int64,size(DT1corner,1)),size(mdbm.posinterp,2),size(mdbm.posinterp,2))
  A+=transpose(A)
  #A=full(A)

  M=(A*A).&A

  x,y=findn(triu(M))
  DT2=Array{Int64}(0,3)
  for k=1:length(x)
    z=find(A[:,x[k]].&A[:,y[k]])
    DT2=vcat(DT2,hcat(x[k] *ones(Int64,size(z)),y[k]*ones(Int64,size(z)),z))
  end

  DT2=mapslices(sort,DT2,2)
  DT2=unique(DT2,1)

  #mdbm.DT1=vcat(DT2[:,[1,2]],DT2[:,[2,3]],DT2[:,[3,1]])
  mdbm.DT2=DT2
  end

end

function mdbm_Gadflylayer(mdbm::mdbm_object,topdim::Int64=1)

if isempty(mdbm.posinterp)
  println("There is no interpolated points!")
  return
end
xsc  = Scale.x_continuous(minvalue=mdbm.ax[1][1], maxvalue=mdbm.ax[1][end])
ysc  = Scale.y_continuous(minvalue=mdbm.ax[2][1], maxvalue=mdbm.ax[2][end])
 if topdim==0 || isempty(mdbm.DT1)
   if length(mdbm.Nax)==1
     x0=mdbm.posinterp[1,:]
     y0=mdbm.posinterp[1,:]*0
   elseif length(mdbm.Nax)>1
     #~~~~~~~~~~~ interploated points ~~~~~~~~~~~~
     x0=mdbm.posinterp[1,:]
     y0=mdbm.posinterp[2,:]
   end
   layer0 = layer(x=x0, y=y0,  Geom.point,style(highlight_width=0.4pt))
   # plot(layer0,xsc,ysc) - not working ???
   return layer0
 elseif !isempty(mdbm.DT1)

   x1=mdbm.posinterp[1,mdbm.DT1[:,1]]
   y1=mdbm.posinterp[2,mdbm.DT1[:,1]]
   x2=mdbm.posinterp[1,mdbm.DT1[:,2]]
   y2=mdbm.posinterp[2,mdbm.DT1[:,2]]

   layer1 = layer(x=x1, y=y1, xend=x2, yend=y2, Geom.segment)
   #plot(layer1,xsc,ysc) - not working ???
   return layer1
 end
end


function mdbm_plot(mdbm::mdbm_object,topdim::Int64=2)
topdim=min(topdim,mdbm.Ndim-mdbm.Ncodim)
if isempty(mdbm.posinterp)
  println("There is no interpolated points!")
  return
end
lims=[(mdbm.ax[k][1], mdbm.ax[k][end]) for k=1:length(mdbm.ax)]

if topdim==0 || isempty(mdbm.DT1)
   if length(mdbm.Nax)==1
     x0=mdbm.posinterp[1,:]
     y0=mdbm.posinterp[1,:]*0
     p= scatter(x0,y0,xlims = lims[1],ylims = (-1,1),markersize =3)
     return p
   elseif length(mdbm.Nax)==2
     #~~~~~~~~~~~ interploated points ~~~~~~~~~~~~
     x0=mdbm.posinterp[1,:]
     y0=mdbm.posinterp[2,:]
     p= scatter(x0,y0,xlims = lims[1],ylims = lims[2],markersize =3)
     return p
   elseif length(mdbm.Nax)>=3
     #~~~~~~~~~~~ interploated points ~~~~~~~~~~~~
     x0=mdbm.posinterp[1,:]
     y0=mdbm.posinterp[2,:]
     z0=mdbm.posinterp[3,:]
     p= scatter(x0,y0,z0,xlims = lims[1],ylims = lims[2],zlims = lims[2],markersize =2)
     return p
   end

elseif topdim==1
  if isempty(mdbm.DT1)
    println("It must be a line object and the line connections must be detected first! Use DTconnect!(mdbm) function")
    return
  end
  #TODO: remove the arrowhead
  if length(mdbm.Nax)==1 #maybe this section is pointles "plotting lines for one-parameter problem
     x1=mdbm.posinterp[1,mdbm.DT1[:,1]]
     y1=x1*0
     x2=mdbm.posinterp[1,mdbm.DT1[:,2]]
     y2=x2*0
     p=quiver(x1, y1, quiver=(x2-x1,y2-y1))
     return p
  elseif length(mdbm.Nax)>=2  #TODO: 3D is not working jet!!!
    x1=mdbm.posinterp[1,mdbm.DT1[:,1]]
    y1=mdbm.posinterp[2,mdbm.DT1[:,1]]
    x2=mdbm.posinterp[1,mdbm.DT1[:,2]]
    y2=mdbm.posinterp[2,mdbm.DT1[:,2]]
    p=quiver(x1, y1, quiver=(x2-x1,y2-y1))
    return p
  elseif length(mdbm.Nax)>=3
    x1=mdbm.posinterp[1,mdbm.DT1[:,1]]
    y1=mdbm.posinterp[2,mdbm.DT1[:,1]]
    z1=mdbm.posinterp[3,mdbm.DT1[:,1]]
    x2=mdbm.posinterp[1,mdbm.DT11[:,2]]
    y2=mdbm.posinterp[2,mdbm.DT1[:,2]]
    z2=mdbm.posinterp[3,mdbm.DT1[:,2]]
    p=quiver(x1, y1, z1, quiver=(x2-x1,y2-y1,z2-z1))
    return p
  end
 elseif topdim==2
   if isempty(mdbm.DT2)
     println("It must be a 'surface' and the toqurface connections must be detected first! Use DTconnect!(mdbm) function")
     return
   end
   println("This part is not ready jet")
 end
end


end
