Ntae_381_eye <- data.frame(ED=4.508, AL=3.1395)
Pmol_753_eye <- data.frame(ED=4.8, AL=3.425)
Pmol_752_eye <- data.frame(ED=5.225, AL=3.9)
eyelist <- list(Ntae_381_eye=Ntae_381_eye ,Pmol_753_eye=Pmol_753_eye ,Pmol_752_eye=Pmol_752_eye)
for (i in 1:length(eyelist)) {
	sph_est <- retinal_arclen(eyelist[[i]]$ED)
	el_est <- semi_ellipse_perimeter(a=eyelist[[i]]$ED, b=eyelist[[i]]$AL)/2
	print(names(eyelist)[i])
	print(paste('ellipse', el_est))
	print(paste('sphere', sph_est))
}
