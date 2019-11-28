import QuantEcon
import Interpolations
import Debugger
import PyPlot
import StatsFuns
import Distributions
import NLsolve
import HDF5

using Distributed
@everywhere using QuantEcon
@everywhere using Interpolations
@everywhere using Debugger
@everywhere using PyPlot
@everywhere using StatsFuns
@everywhere using Distributions
@everywhere using NLsolve
@everywhere using HDF5



type parameters_type
	M::Float64
	ddelta::Float64
	aalpha_k::Float64
	aalpha_l::Float64
	ggamma::Float64
	B::Float64
	w::Float64
	ttheta::Float64
	f::Float64
	C::Float64

	rrho::Float64
	ssigma::Float64
	Nz::Int64
	Pz::Array{Float64,2}
	z_grid::Array{Float64}

	Nk::Int64
	K_grid::Array{Float64}
	Nk_prime::Int64
	K_prime_grid::Array{Float64}
	new_points::Int64

	prec_VFI::Float64
	prec_sol::Float64
	prec_sol_B::Float64
	accelerator_flag::Int64
	accelerator_step::Int64
	maxiter::Int64

	K_init::Float64
	T_simul::Int64
	B_target::Float64
end

function define_params(w)
	M=0.95
	ddelta=0.1

	aalpha_k=0.3
	aalpha_l=0.65
	ggamma=aalpha_k/(1-aalpha_l)
	f=0.011

	B=(1-aalpha_l)*(w/aalpha_l)^(aalpha_l/(aalpha_l-1))

	K_ss=0.75
	suppl_const=(1-aalpha_l)*(w/aalpha_l)^(aalpha_l/(aalpha_l-1))
	C=(1/M-1+ddelta)/(suppl_const*ggamma*K_ss^(ggamma-1)) #normalization constant (to get K_ss=0.75)


	Nz=15
	rrho=0.95
	ssigma=0.02
	N_std=3
	bar_z=0.0
	a_t=tauchen(Nz, rrho, ssigma,bar_z,N_std) #Standard QuantEcon function
	Pz=a_t.p  #transition matrix
	z_grid=exp(a_t.state_values) #grid for TFP

	Nk=401
	K_min=0.001
	K_max=20.0
	#K_grid=linspace(K_min,K_max,Nk)
	K_grid=Array(Float64,Nk)
	step=(K_max/K_min)^(1/(Nk-1))
	K_grid[1]=K_min
	for iK=2:Nk
		K_grid[iK]=K_grid[iK-1]*step
	end

	new_points=9 #how many points to throw in each interval between neighbouring grid nodes of S_grid
	Nk_prime=(1+new_points)*(Nk-1)+1

	K_prime_grid=Array(Float64,Nk_prime)
	l=1
	for iK=1:Nk-1
		K_prime_grid[l:l+new_points+1]=linspace(K_grid[iK],K_grid[iK+1],new_points+2)
		l=l+new_points+1
	end


	#adj costs
	ttheta=0.5





	prec_VFI=0.000001
	prec_sol=0.000000001
	prec_sol_B=0.1
	accelerator_flag=1
	accelerator_step=20
	maxiter=1000

	T_simul=5000000
	K_init=K_grid[1]
	B_target=1.0

	p=parameters_type(M,ddelta,aalpha_k,aalpha_l,ggamma,B,w,ttheta,f,C,
						rrho,ssigma,Nz,Pz,z_grid,
						Nk,K_grid,Nk_prime,K_prime_grid,new_points,
						prec_VFI,prec_sol,prec_sol_B,accelerator_flag,accelerator_step,maxiter,
						K_init,T_simul,B_target)
	return p
end



function profit_func(K::Float64,z::Float64,p::parameters_type)
	profit=p.C*K^p.ggamma*p.B*z^(1/(1-p.aalpha_l))-p.f
	return profit
end

function sales_func(K::Float64,z::Float64,p::parameters_type)
	sales=p.C*K^p.ggamma*p.B*z^(1/(1-p.aalpha_l))/(1-p.aalpha_l)
	return sales
end


function VFI(V_init::Array{Float64},p::parameters_type)

	#Optimal policies
	K_prime_arr=Array(Float64,p.Nk,p.Nz)
	iK_prime_arr=Array(Int64,p.Nk,p.Nz)
	Div_arr=Array(Float64,p.Nk,p.Nz)
	Exit_arr=zeros(Int64,p.Nk,p.Nz)

	V=V_init


	diff=1 #stopping criterion: maximum(abs(V_previos_step-V_current_step))<prec_VFI
	iter=1 #initiate iterator

	#main loop
	println("Value function iteration...")
	while diff>p.prec_VFI && iter<p.maxiter

		acceleration=0
		if mod(iter,p.accelerator_step)!=1 && p.accelerator_flag==1 && iter>2
			acceleration=1
		end

		V_prev=copy(V)
		EV=V_prev*p.Pz'

		EV_i=interpolate((p.K_grid,p.z_grid),EV,Gridded(Linear())) #linear interpolation of EV to find optimal policies more precisely
		EV_interp=EV_i[p.K_prime_grid,p.z_grid]

		if acceleration==0 #regular step of VFI

			for iK=1:p.Nk

				K=p.K_grid[iK]
				for iz=1:p.Nz
					z=p.z_grid[iz]

					profit=profit_func(K,z,p)
					inv=p.K_prime_grid-(1-p.ddelta)*K
					adj=p.ttheta*(p.K_prime_grid/K-1).^2*K
					div=profit-inv-adj-p.f
					V_cur=div+p.M*max(EV_interp[:,iz],0)

					pos=indmax(V_cur)
					V[iK,iz]=V_cur[pos]
					K_prime_arr[iK,iz]=p.K_prime_grid[pos]
					Div_arr[iK,iz]=div[pos]
					iK_prime_arr[iK,iz]=pos
					if max(EV_interp[pos,iz],0)==0
						Exit_arr[iK,iz]=1
					else
						Exit_arr[iK,iz]=0
					end

				end
			end




		else
			for iK=1:p.Nk
				K=p.K_grid[iK]
				for iz=1:p.Nz
					z=p.z_grid[iz]
					V[iK,iz]=Div_arr[iK,iz]+p.M*max(EV_i[K_prime_arr[iK,iz],z],0)
				end
			end
		end

		diff=maximum(abs(V-V_prev))

		if acceleration==1
			#println("Iteration:",iter,"  Diff: ", round(diff/p.prec_VFI-1,2))
			diff=1
		else
			println("Iteration:",iter,"  Diff: ", round(diff/p.prec_VFI-1,2))
		end
		iter+=1
	end


println("Done!")

return V,iK_prime_arr,K_prime_arr,Div_arr,Exit_arr

end


function run_VFI(read_initial_guess_flag,w)
	p=define_params(w)


	if read_initial_guess_flag==1
		V_init=h5read("results_VFI.h5","V") #read from file and then interpolate to get a good initial guess
		K_grid_init=h5read("results_VFI.h5","K_grid")
		z_grid_init=h5read("results_VFI.h5","z_grid")

		V_init_i=interpolate((K_grid_init,z_grid_init),V_init,Gridded(Linear()))

		V_init=V_init_i[p.K_grid,p.z_grid]
	else
		V_init=zeros(p.Nk,p.Nz)
	end

	V,iK_prime,K_prime,Div,Exit=VFI(V_init,p)


	h5open("results_VfI.h5","w") do file
		   write(file, "V", V)
		   write(file, "iK_prime", iK_prime)
		   write(file, "K_prime", K_prime)
		   write(file, "Div", Div)
		   write(file, "Exit", Exit)
		   write(file, "K_grid", convert(Array{Float64},p.K_grid))
		   write(file, "z_grid", convert(Array{Float64},p.z_grid))
	end

	return V,iK_prime,K_prime,Div,Exit,p
end



function lr_prob(p::parameters_type)
	Pz_lr=p.Pz^1000
	long_run_prob=Pz_lr[1,:]
	return long_run_prob
end


function full_transition(iK_prime::Array{Int64},p::parameters_type)
	T=zeros(Float64,p.Nk*p.Nz,p.Nk*p.Nz)

	l=1
	for iz=1:p.Nz
		for iK=1:p.Nk
			pos=1
			for iz_prime=1:p.Nz

					t=floor(Int64,iK_prime[iK,iz]/(p.new_points+1))

					pos_down=t+1
					pos_up=t+2

					prob_down=(p.K_prime_grid[iK_prime[iK,iz]]-p.K_grid[pos_down])/(p.K_grid[pos_up]-p.K_grid[pos_down])
					prob_up=(p.K_grid[pos_up]-p.K_prime_grid[iK_prime[iK,iz]])/(p.K_grid[pos_up]-p.K_grid[pos_down])

					T[pos_down+(iz_prime-1)*p.Nk,l]=p.Pz[iz,iz_prime]*prob_down
					T[pos_up+(iz_prime-1)*p.Nk,l]=p.Pz[iz,iz_prime]*prob_up


			end
			l+=1
		end
	end

	return T
end


#@debug
function dstrb_no_entry(iK_prime::Array{Int64},p::parameters_type)

	T=full_transition(iK_prime,p)
	a=T^150
	mu_stat=a[:,1]

	dens=reshape(mu_stat,p.Nk,p.Nz)/sum(mu_stat)

	fig = figure("Stationary distribution of firms",figsize=(10,10))
	ax = fig[:add_subplot](1,1,1, projection = "3d")
	ax[:plot_surface](repmat(p.K_grid,1,p.Nz), repmat(p.z_grid',p.Nk,1), dens, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
	#ax[:set_ylim]([0,5])
	ylabel("Shock")
	xlabel("Capital")
	title("Stationary distribution of firms")

	return mu_stat,dens,T
end

@debug function dstrb_with_entry(iK_prime::Array{Int64},Exit::Array{Int64},p::parameters_type,w::Float64)

	T=full_transition(iK_prime,p)
	l=1
	for iz=1:p.Nz
		for iK=1:p.Nk
			if Exit[iK,iz]==1
				T[:,l]=0]
			end
			l+=1
		end
	end

	suppl=p.Pz^250
	lr_prob=suppl[1,:]
	suppl=0
	ggamma=zeros(p.Nk*p.Nz)
	for iz=1:p.Nz
		ggamma[(iz-1)*p.Nk+1]=lr_prob[iz]
	end


	mu_stat=(eye(p.Nk*p.Nz)-T)\ggamma
	E=1/sum(mu_stat)
	mu_stat=mu_stat*E

	labor=Array(Float64,p.Nk,p.Nz)
	profit=Array(Float64,p.Nk,p.Nz)
	for iK=1:p.Nk
		for iz=1:p.Nz
			labor[iK,iz]=labor_demand_func(w,p.K_grid[iK],p.z_grid[iz],p)
			profit[iK,iz]=profit_func(p.K_grid[iK],p.z_grid[iz],p)
		end
	end
	dens=reshape(mu_stat,p.Nk,p.Nz)/sum(mu_stat)

	l_demand=sum(dens.*labor)
	profit_av=sum(profit.*dens)
	@bp
	B=l_demand/w^0.5

	entry=sum(E*ggamma)
	entry_cost=sum(lr_prob.*V[1,:])


	fig = figure("Stationary distribution of firms",figsize=(10,10))
	ax = fig[:add_subplot](1,1,1, projection = "3d")
	ax[:plot_surface](repmat(p.K_grid,1,p.Nz), repmat(p.z_grid',p.Nk,1), dens, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
	ylabel("Shock")
	xlabel("Capital")
	title("Stationary distribution of firms")

	return mu_stat,dens,T,E,entry,entry_cost,l_demand,B,profit_av
end


function labor_demand_func(w::Float64,k::Float64,z::Float64,p::parameters_type)
	return (w/(p.C*z*k^p.aalpha_k))^(1/(p.aalpha_l-1))
end
