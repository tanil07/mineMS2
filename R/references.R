
###MACRO REFERENCES

HIGH_MASS_SYMBOL <- NA_character_
HIGH_MASS_INDICATOR <- NA_character_
HIGH_MASS_LEGEND <- "high m/z"
UNKNOWN_FORMULA_INDICATOR <- NA_real_
UNKNOWN_FORMULA_LEGEND <- "unknown formula"


#################
#ATOMS REFERENCE#
#################

###Reference for atoms
tabAtoms <- function(){
	n_atoms <- c("H", "C", "O", "N", "F", "Cl", "Br", "I", "Si",
				 "S", "P", "Se")
	m_atoms <- c(1.00783, 12, 15.994915, 14.003074, 18.998403,
				 34.0968853, 78.918336, 126.904477, 27.976928, 31.972071,
				 30.973763, 73.922477)
	val_atoms <- c(1, 4, 2, 3, 1, 1, 1, 1, 4, 2, 3, 6)
	halogens <- c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)
	return(data.frame(name = n_atoms, mass = m_atoms, valence = val_atoms, halogen = halogens,
					  stringsAsFactors = FALSE))
}

###Construct a reference table form a lsit of atoms.
makeDfValidAtoms <- function (clist, mz = NULL)
{
	tAtoms <- tabAtoms()
	ref <- names(clist)
	pok <- ref %in% tAtoms$name
	if (any(!pok)) {
		stop("invalid atoms furnished :", paste(ref[!pok], collapse = ",",
												sep = " "))
	}
	pRef <- match(ref[pok], tAtoms$name)
	nmaxAtoms <- ceiling(mz/tAtoms$mass[pRef])
	vnum <- NULL
	if (class(clist) == "list") {
		vnum <- unlist(clist)
	}
	else {
		vnum <- clist
	}
	dfres <- data.frame(name = ref, mass = tAtoms$mass[pRef],
						valence = tAtoms$valence[pRef])
	if (!is.null(mz))
		dfres$nums <- apply(rbind(unlist(clist), nmaxAtoms),
							2, min)
	return(dfres)
}


####################
#ADDUCTS REFERENCES#
####################

tabAdducts<-function(charge=NULL){
	tAdd <- data.frame(name=c('','H','Na','K','NH4',
							  '2Na-H','2K-H','2NH4-H',
							  '2H','2Na','2K','2NH4'),
					   charge=c(1,1,1,1,1,1,1,1,2,2,2,2),
					   quasi=c(0,1,1,1,1,0,0,0,0,0,0,0),
					   massdiff=c(0,1.007276,22.989218,38.963158,18.033823,
					   		   44.97116,76.91904,35.06037,
					   		   2.014552,45.97844,77.92632,36.06765),
					   stringsAsFactors = FALSE
	)
	# ###TO DEBUG ONLY
	# tAdd <- tAdd[-c(5,8,12),]
	if(!is.null(charge)){
		tAdd <- tAdd[which(tAdd$charge<=charge),]
		if(charge<0){
			tAdd$charge <- (-tAdd$charge)
			tAdd$massdiff <- (-tAdd$massdiff)
		}
	}
	tAdd
}

makeDfValidAdducts<-function(adducts,charge=2){

	tAdducts <- tabAdducts(abs(charge))
	pok <- adducts %in% tAdducts$name
	if(any(!pok)){
		stop("invalid adduct furnished :",paste(adducts[!pok],collapse = ",",sep=" "))
	}
	pRef <- match(adducts[pok],tAdducts$name)
	dfres <- data.frame(name=adducts,massdiff=tAdducts$massdiff[pRef],
						charge=tAdducts$charge[pRef],stringsAsFactors = FALSE)
	return(dfres)
}

####################
#LOSSES REFERENCES#
####################


#' Default table of penalized loss
#'
#' The losses are penalized by the scoring function,
#' they are extracted form the paper by Bocker et al 2010. This includes both very
#' frequent losses, which will lead ot irllevant patterns and forbidden losses.
#'
#' @return A data frame with 2 columns the names of the losses and the
#' formula of the losses. Only the second column is used the first column is informative.
#' @export
#'
#' @examples
#' penalizedLossesDefault()
penalizedLossesDefault <- function() {
	structure(
		list(
			name = c(
				"Water",
				"Methane",
				"Ethene",
				"Ethine",
				"Butene",
				"Pentene",
				"Benzene",
				"Formaldehyde",
				"Methanol",
				"Carbon Monoxide",
				"Formic acid",
				"Carbon dioxide",
				"Acetic acid",
				"Ketene",
				"Propionic acid",
				"Malonic acid",
				"Malonic anhydride",
				"Pentose equivalent",
				"Deoxyhexose equivalent",
				"Hexose equivalent",
				"Hexuronic equivalent acid",
				"Ammonia",
				"Methylamine",
				"Methylimine",
				"Trimethylamine",
				"Cyanic Acid",
				"Urea",
				"Phosphonic acid",
				"Phosphoric acid",
				"Metaphosphoric acid",
				"Dihydrogen vinyl phosphate",
				"Hydrogen sulfide",
				"Sulfur",
				"Sulfur dioxide",
				"Sulfur trioxide",
				"Sulfuric acide",
				"Dicarbon monoxide",
				"Tetracarbon monoxide",
				"Unsaturated cyclopropane",
				"Unsaturated cyclopentane",
				"Unsaturated cycloheptane"
			),
			formula = c(
				"H2O",
				"CH4",
				"C2H4",
				"C2H2",
				"C4H8",
				"C5H8",
				"C6H6",
				"CH2O",
				"CH4O",
				"CO",
				"CH2O2",
				"CO2",
				"C2H4O2",
				"C2H2O",
				"C3H6O2",
				"C3H4O4",
				"C3H2O3",
				"C5H8O4",
				"C6H10O4",
				"C6H10O5",
				"C6H8O6",
				"NH3",
				"CH5N",
				"CH3N",
				"C3H9N",
				"CHNO",
				"CH4N2O",
				"H3PO3",
				"H3PO4",
				"HPO3",
				"C2H5O4P",
				"H2S",
				"S",
				"SO2",
				"SO3",
				"H2SO4",
				"C2O",
				"C4O",
				"C3H2",
				"C5H2",
				"C7H2"
			)
		),
		.Names = c("name", "formula"),
		class = "data.frame",
		row.names = c(NA,-36L)
	)
}

