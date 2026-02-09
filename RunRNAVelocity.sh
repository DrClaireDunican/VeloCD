echo "Hello, thank you for using VeloCD. It was developed by Dr Claire Dunican as a part of her PhD and post-doctoral research. Before we start I will ask you a few questions. 
What is the full location of your expression files? Please type it in the format: /user/Xfolder/Yfolder"

read MainExpressionFileLoc

until [ -d "$MainExpressionFileLoc" ]
do
    echo "That is not a valid location. What is the full location of your expression files? Please type it in the format: /user/Xfolder/Yfolder"
    read MainExpressionFileLoc
done

echo "Where is the source code of VeloCD stored on your computer? Please type it in the format: /user/Xfolder/Yfolder"

read CodeLocation

until [ -d "$CodeLocation" ]
do
    echo "That is not a valid location. What is the full location of the VeloCD source code? Please type it in the format: /user/Xfolder/Yfolder"
    read CodeLocation
done

echo "Please select the mode you would like use. Mode 1 (type: M1): for all the sets of files in a directory. Mode 2 (type: M2) for one set of files in the directory (specified above). Each set consists of one or two metadata files (Metadata_*.csv), a spliced transcript expression file (Spliced_*.csv) and unspliced transcript expression file (Unspliced_*.csv), where * is a common substring per set. For example, Spliced_SixteenGeneSignature.csv, where *=SixteenGeneSignature. If you have multiple metadata files per set, please ensure that the second file is titled in the format: MetadataAlt_*.csv."

read RunMode

until [[ "$RunMode" = @("M1"|"M2") ]]
do
    echo "That was not a valid mode, please type M1 or M2."
    read RunMode
done

if [ "$RunMode" = "M2" ]
then
    echo "What is the common substring of your set?"
    read CommonSubstring
fi

echo "PCA will be run as a standard part of VeloCD."

echo "How many PCA components would you like to generate. Please type any integer greater than 1. If you only have two genes in your dataset then please type 2."

read PCANum

until [[ "$PCANum" =~ ^[0-9]+$ ]] && [[ "$PCANum" -gt 1 ]];
do
    echo "That was not a valid answer: either non-numerical or less than the minimum value, pleas try again. How many PCA components would you like to generate. Please type any integer greater than 1."
    read PCANum
done

if [ "$PCANum" = "2" ]
then
    PCAN="2"
elif [ "$PCANum" = "3" ]
then
    PCAN="3"
else
    echo "How many PCA components would you like to plot? PLease type 2 or 3."
    read PCANComp
    until [[ $PCANComp = @(2|3) ]]
    do
        echo "That was not 2 or 3. How many PCA components would you like to plot? PLease type 2 or 3."
        read PCANComp
    done
    if [ "$PCANComp" = "2" ]
    then
        PCAN="2"
    else
        PCAN="3"
    fi
fi

echo "Would you also like to embed your fate map arrows using PCA? Please type yes or no."

read PCAYesNo

until [[ "$PCAYesNo" = @("yes"|"no") ]];
do
    echo "Invalid answer. Would you also like to embed your fate map arrows using PCA? Please type yes or no."
    read PCAYesNo
done

if [ "$PCAYesNo" = "yes" ]
then
    if [[ $PCANum -gt 2 ]] && [[ $PCAN == 3 ]]
    then
        echo "Would you like to have 3-dimensional PCA-based fate maps? Please type yes or no."
        read ThreeDPCA
        until [[ "$ThreeDPCA" = @("yes"|"no") ]];
        do
            echo "Invalid answer. Would you like to have 3-dimensional PCA-based fate maps? Please type yes or no."
            read ThreeDPCA
        done
    else
        ThreeDPCA="no"
    fi
    echo "Would you like the neighbourhood search to weight the dimensions by the percentage of expression variation explained by each Principal Component? Please type yes or no."
    read AdjustYesNo
    until [[ "$AdjustYesNo" = @("yes"|"no") ]];
    do
        echo "Invalid answer. Please type yes or no."
        read AdjustYesNo
    done
    if [ "$AdjustYesNo" = "yes" ]
    then
        AdjustPCANDS="True"
    else 
        AdjustPCANDS="False"
    fi
else
    ThreeDPCA="no"
fi


echo "Would you like to embed your fate maps using tSNE? Please type yes or no. 
Please note that tSNE is one of the slower embedding algorithms and may limit the speed of this process."

read tSNEYesNo

until [[ "$tSNEYesNo" = @("yes"|"no") ]];
do
    echo "Invalid answer. Would you like to embed your fate maps using tSNE? Please type yes or no."
    read tSNEYesNo
done

echo "Would you like to embed your fate maps in UMAP? Please type yes or no."

read UMAPYesNo

until [[ "$UMAPYesNo" = @("yes"|"no") ]];
do
    echo "Invalid answer. Would you like to embed your fate maps using UMAP? Please type yes or no."
    read UMAPYesNo
done


if [ "$UMAPYesNo" = "yes" ]
then
    echo "What distance metric would you like to use for this embedding? Options include euclidean, manhattan, chebyshev and minkowski. For more information, please see: the metric section of https://umap-learn.readthedocs.io/en/latest/parameters.html"
    read UMAPDist
    until [[ "$UMAPDist" = @("euclidean"|"manhattan"|"chebyshev"|"minkowski") ]];
    do
        echo "Invalid answer. What distance metric would you like to use for this embedding? The options are: euclidean, manhattan, chebyshev and minkowski. For more information, please see: the metric section of https://umap-learn.readthedocs.io/en/latest/parameters.html"
        read UMAPDist
    done
    echo "What minimum distance would you like to use for this embedding method? This gives the minimum distance between points in the embedding. Values range between 0.1 and 0.99. The default is 0.1. For more information, please see: the min_dist section of https://umap-learn.readthedocs.io/en/latest/parameters.html"
    read UMAPDistNum
    echo "Would you also like 3-dimensional UMAP embeddings? Please type yes or no."
    read UMAP3DYesNo
    until [[ "$UMAP3DYesNo" = @("yes"|"no") ]];
    do
        echo "Invalid answer. Would you also like 3-dimensional UMAP embeddings? Please type yes or no."
        read UMAP3DYesNo
    done
else
    UMAPDist="euclidean"
    UMAPDistNum="0.1"
    UMAP3DYesNo="no"
fi


if [ "$UMAPYesNo" = "yes" ] || [ "$tSNEYesNo" = "yes" ]
then
    echo "Would you like to embed your fate maps across multiple embedding number of neighbour hyperparameter values? Please type yes or no."
    read RunMultiEmbed
    until [[ "$RunMultiEmbed" = @("yes"|"no") ]]
    do
        echo "That was not a valid response, please type yes or no."
        read RunMultiEmbed
    done
    if [ $RunMultiEmbed = "yes" ]
    then
        echo "What range of embedding number of neighbours would you like to use? Please give the minimum value, this has to be at least 2."
        read Embed1
        until [[ "$Embed1" =~ ^[0-9]+$ ]] && [[ "$Embed1" -gt 1 ]];
        do
            echo "That was not a valid response, either non-numerical or less than the minimal value. What range of embedding number of neighbours would you like to use? Please give the minimum value, this has to be at least 2."
            read Embed1
        done
        echo "Please give a maximum value, this has to be maxmimum the number of samples minus 1."
        read Embed2
        until [[ "$Embed2" =~ ^[0-9]+$ ]] && [[ "$Embed2" -gt "$Embed1" ]]; #issue allows 7 if prev is 5 but not 12 or 37
        do
            echo "That was not a valid response: either non-numerical or less than the minimal value. Please give a maximum value, this has to be maxmimum the number of samples minus 1."
            read Embed2
        done
    else
        echo "Which embedding value would you like to use?"
        read Embed1
        until [[ "$Embed1" =~ ^[0-9]+$ ]] && [[ "$Embed1" -gt 1 ]];
        do
            echo "That was not a valid response, either non-numerical or less than the minimal value. Please type a number."
            read Embed1
        done
        Embed2="10"
    fi
else
    if [ "$UMAPYesNo" = "no" ] || [ "$tSNEYesNo" = "no" ] || [ "$PCAYesNo" = "yes" ]
    then
        Embed1="NotApplicable"
        Embed2="NotApplicable"
        RunMultiEmbed="no"
   else
        Embed1=2
        Embed2=3
        RunMultiEmbed="no"
   fi
fi

if [ "$UMAPYesNo" = "yes" ] || [ "$tSNEYesNo" = "yes" ] || [ "$PCAYesNo" = "yes" ]
then
    echo "It is highly recommended that VeloCD is run with a range of transition probability number of neighbour values. What range of transition probabilities would you like to use? Please give the minimum value, this has to be at least 2."
    read TPMin
    until [[ "$TPMin" =~ ^[0-9]+$ ]] && [[ "$TPMin" -gt 1 ]];
    do
        echo "Invalid response: either less than 2 or not a number. What range of transition probabilities would you like to use? Please give the minimum value, this has to be at least 2."
        read TPMin
    done
    echo "Please give the maximum value plus 1, this has to be less than the number of samples."
    read TPMax
    until [[ "$TPMax" =~ ^[0-9]+$ ]] && [[ "$TPMax" -gt "$TPMin" ]]; #issue allows 7 if prev is 5 but not 12 or 37
    do
        echo "That was not a valid response: either non-numerical or less than the minimal value.  Please give the maximum value plus 1, this has to be less than the number of samples."
        read TPMax
    done
    echo "Please give the step between these two values. i.e. 1 for a range 13 to 16 would be 13, 14, 15 and 16. 2 would give 13 and 15."
    read TPStep
    until [[ "$TPStep" =~ ^[0-9]+$ ]] && [[ "$TPStep" -gt 0 ]];
    do
        echo "That was not a valid response: either non-numerical or less than the minimal value. Please give the step between these two values. i.e. 1 for a range 13 to 16 would be 13, 14, 15 and 16. 2 would give 13 and 15."
        read TPStep
    done
    echo "What minimum probability threshold would you like to use for prediction? Please type a value between 0 and 1. The default is 0.5, which is preferable when there are two groups."
    read ProbabilityThreshold
else 
    TPMin=2
    TPMax=4
    TPStep=1
    ProbabilityThreshold="0.5"
fi

if [ "$UMAPYesNo" = "yes" ] || [ "$tSNEYesNo" = "yes" ] || [ "$PCAYesNo" = "yes" ] 
then
    if [ "$UMAP3DYesNo" = "yes" ] || [ "$ThreeDPCA" = "yes" ]
    then
        if [[ $PCANum -gt 2 ]] && [[ $PCAN == 3 ]]
        then
            echo "Would you like to make GIFs of the 3-dimensional plots (PCA and/or UMAP)? Please type yes or no. This process will increase the run-time of this method."
            read GIFY
            until [[ "$GIFY" = @("yes"|"no") ]];
            do
                echo "Invalid Response. Would you like to make GIFs of the 3-dimensional plots (PCA and/or UMAP)? Please type yes or no. This process will increase the run-time of this method."
                read GIFY
            done
            if [ "$GIFY" = "yes" ]
            then
                GIFRUN="True"
            else
                GIFRUN="False"
            fi
        else
            GIFRUN="False"
        fi
    else 
        GIFRUN="False"
    fi
else
    GIFRUN="False"
fi


echo "Would you like to shape the 'sample' points on the plots by your metadata variable? Please type yes or no. This is not recommended when you have more than 3 groups. If only one metadata file is provided, these will match the colours given. If a a second is provided, these markers will be determined using the 'Group' column in the second file."

read ShapeYesNo

until [[ "$ShapeYesNo" = @("yes"|"no") ]];
do
    echo "Invalid response. Would you like to shape the points by your metadata variable? Please type yes or no. This is not recommended when you have more than 3 groups."
    read ShapeYesNo
done

if [ "$ShapeYesNo" = "yes" ]
then
    ShapePoint="True"
else
    ShapePoint="False"
fi

echo "Which colourmap would you like to use to colour the sample points on the plots? Please type in the name of one from matplotlib: https://matplotlib.org/stable/tutorials/colors/colormaps.html. I personally recommend Set1."

read ColorMap


echo "What type of feature is your meta-feature of interest? If they are categorical plese type: strings. If they are numerical please type: integers. For non-numerical features please ensure your names don't have special characters."

read Typeof1

until [[ "$Typeof1" = @("strings"|"integers") ]];
do
    echo "Invalid reponse. What type of feature is your meta-feature of interest? If they are categorical plese type: strings. If they are numerical please type: integers"
    read Typeof1
done

echo "Do you have two unique sample class labels in your primary metadata file (over which you are predicting)? Please type yes or no. If you have more than two then type no."

read NumGroups

until [[ "$NumGroups" = @("yes"|"no") ]];
do
    echo "Invalid response. Please type yes or no."
    read NumGroups
done 

if [ "$NumGroups" = "yes" ]
then
    if [ "$UMAPYesNo" = "yes" ] || [ "$tSNEYesNo" = "yes" ] || [ "$PCAYesNo" = "yes" ]
    then
        echo "Would you like sensitivity and specificity values to be calculated for all runs using all the samples in each set? Please type yes or no"
        read SensSpec
        until [[ "$SensSpec" = @("yes"|"no") ]];
        do
            echo "Invalid response. Please type yes or no."
            read SensSpec
        done
        if [ "$SensSpec" = "yes" ]
        then 
            echo "What is the feature label that you would like to use to calculate sensitivity? These features are present in your primary metadata file."
            read Feature1
            echo "What is the feature label that you would like to use to calculate specifity?"
            read Feature2
        else
            Feature1="Fatality"
            Feature2="Survivor"
        fi
    else 
        SensSpec="no"
        Feature1="Fatality"
        Feature2="Survivor"
    fi
else 
    SensSpec="no"
    Feature1="Fatality"
    Feature2="Survivor"
fi

echo "Do you have multiple metadata files per run? This second file would contain an addiitonal variable for fate map visualisation (groups that define the shape of sample points) and to predict upon (optional, you will be asked if you want to do this in a later question). PLlease type yes or no."

read SeconMe

until [[ "$SeconMe" = @("yes"|"no") ]];
do
    echo "Do you have multiple metadata files per run? This second file would contain an addiitonal variable for fate map visualisation (groups that define the shape of sample points) and to predict upon (optional, you will be asked if you want to do this in a later question). PLlease type yes or no."
    read SeconMe
done

if [ "$SeconMe" = "yes" ]
then
    SecondMeta="True"
    echo "What type of feature is your second meta-feature of interest? If they are categorical plese type: strings. If they are numerical please type: integers"
    read Typeof2
    until [[ "$Typeof2" = @("strings"|"integers") ]];
    do
        echo "Invalid reponse. What type of feature is your meta-feature of interest? If they are categorical plese type: strings. If they are numerical please type: integers"
        read Typeof2
    done
    if [ "$UMAPYesNo" = "yes" ] || [ "$tSNEYesNo" = "yes" ] || [ "$PCAYesNo" = "yes" ]
    then
        echo "Would like to use this second meta-feature for prediction? Please type yes or no. If you type no this feature will only be used to shape the points on the fate maps"
        read PredSecond
        until [[ "$PredSecond" = @("yes"|"no") ]];
        do
            echo "Invalid reponse. Please type yes or no."
            read PredSecond
        done
        if [ "$PredSecond" = "yes" ]
        then
            echo "Do you have two unique sample class labels in this additional metadata file (over which you are predicting)? Please type yes or no. Type no if you have more than two."
            read NumGroups2
            until [[ "$NumGroups2" = @("yes"|"no") ]];
            do
                echo "Invalid response. Please type yes or no."
                read NumGroups2
            done 
            if [ "$NumGroups2" = "yes" ]
            then
                echo "Would you like sensitivity and specificity values to be calculated for all runs using all the samples in each set? Please type yes or no"
                read SensSpec2
                until [[ "$SensSpec2" = @("yes"|"no") ]];
                do
                    echo "Invalid response. Please type yes or no."
                    read SensSpec2
                done
                if [ "$SensSpec2" = "yes" ]
                then 
                    echo "What is the feature label that you would like to use to calculate sensitivity? These features are present in your primary metadata file."
                    read Feature21
                    echo "What is the feature label that you would like to use to calculate specifity?"
                    read Feature22
                else
                    Feature21="Fatality"
                    Feature22="Survivor"
                fi
            fi
        fi
    else
        Feature21="Fatality"
        Feature22="Survivor"
        SecondMeta="False"
        Typeof2="strings"
        PredSecond="no"    
    fi
else
    Feature21="Fatality"
    Feature22="Survivor"
    SecondMeta="False"
    Typeof2="strings"
    PredSecond="no"
fi

if [ "$Typeof1" = "strings" ]
then
    Typeof1="Strings"
fi

if [ "$Typeof2" = "strings" ]
then
    Typeof2="Strings"
fi


#echo "Would you like to use any additional features?"

#read AditionalFeat

#until [[ "$AditionalFeat" = @("yes"|"no") ]];
#do
#    echo "Invalid response. Would you to use any additional features?"
#    read AditionalFeat
#done

#if [ "$AditionalFeat" = "yes" ]
#then
#    echo "Would you like confidence intervals to be calculated? Please type yes or no."
#    read ConfidenceInterval
#    until [[ "$ConfidenceInterval" = @("yes"|"no") ]];
#    do
#        echo "Invalid response. Would you like confidence intervals to be calculated? Please type yes or no."
#        read ConfidenceInterval
#    done    
#    if [ "$ConfidenceInterval" = "yes" ]
#    then
#        CI="True"
#        echo "How many iterations would you like to do for this calculation? Please type a number."
#        read CIIteration
#        until [[ "$CIIteration" =~ ^[0-9]+$ ]] && [[ $CIIteration -gt 0 ]];
#        do
#            echo "That was not a valid response: either non-numerical or less than the minimal value. How many iterations would you like to do for this calculation? Please type a number."
#            read CIIteration
#        done
#    else
#        CI="False"
#        CIIteration="1000"
#    fi
#    echo "Would you like to perform feature selection based on the contribution of each gene to each Principal Component? Please type yes or no."
#    read PCVAR
#    until [[ "$PCVAR" = @("yes"|"no") ]];
#    do
#        echo "Invalid response. Would you like to perform feature selection based on the contribution of each gene to each Principal Component? Please type yes or no."
#        read PCVAR
#    done 
#    if [ "$PCVAR" = "yes" ]
#    then
#        PCAFS="PCAVar"
#    else
#        PCAFS="NoPCAVar"
#    fi
#    echo "Would you like to perform correlation coefficient-based feature selection? Please type yes or no. This is only to be used when the meta-feature is numerical."
#    read CCFS
#    until [[ "$CCFS" = @("yes"|"no") ]];
#    do
#        echo "Invalid response. Would you like to perform correlation coefficient-based feature selection? Please type yes or no. This is only to be used when the meta-feature is numerical."
#        read CCFS
#    done 
#    if [ "$CCFS" = "yes" ]
#    then
#        CCFSOption="CoorelationCoeff"
#    else
#        CCFSOption="NoCoorelationCoeff"
#    fi
#    echo "Would you like to perform random sampling of the samples? Please type yes or no. This is only recommended when there are hundreds of samples."
#    read RandomSamp
#    until [[ "$RandomSamp" = @("yes"|"no") ]];
#    do
#        echo "Would you like to perform random sampling of the samples? Please type yes or no. This is only recommended when there are hundreds of samples."
#        read RandomSamp
#    done
#    if [ "$RandomSamp" = "yes" ]
#    then
#        RS="True"
#    else
#        RS="False"
#    fi
#else
RandomSamp="no"
RS="False"
ConfidenceInterval="no"
CI="False"
CIIteration="10"
PCVAR="no"
PCAFS="NoPCAVar"
CCFS="no"
CCFSOption="NoCoorelationCoeff"
#fi

wd="${MainExpressionFileLoc}/"

#sub-directories for output

echo "----------------------------------------------------------"

echo "Setting up your output file directories..."

mkdir ${MainExpressionFileLoc}/Results

Main="${MainExpressionFileLoc}/Results"

mkdir $Main/Fate_Maps/
mkdir $Main/Fate_Maps/PCA/
mkdir $Main/Files/
mkdir $Main/Files/Velocity_Values/
mkdir $Main/Files/Predicted_Unspiced/
mkdir $Main/Files/Future_Spliced/

if [ "$PCAYesNo" = "yes" ]
then
    mkdir $Main/Files/PCA/
    mkdir $Main/Files/PCA/Delta_Embedding/
    mkdir $Main/Files/PCA/Self_Correlation_Coefficients/
    mkdir $Main/Files/PCA/Unitary_Vector_X/
    mkdir $Main/Files/PCA/Unitary_Vector_Y/
    mkdir $Main/Files/PCA/Transition_Probability/
    mkdir $Main/Files/PCA/Transition_Probability/Sample_Level/
    mkdir $Main/Files/PCA/Transition_Probability/Group_Level/
    if [ "$ThreeDPCA" = "yes" ]
    then
        mkdir $Main/Files/PCA/Unitary_Vector_Z/
    fi
fi


if [ "$tSNEYesNo" = "yes" ]
then
    mkdir $Main/Fate_Maps/tSNE/
    mkdir $Main/Files/tSNE/
    mkdir $Main/Files/tSNE/Delta_Embedding/
    mkdir $Main/Files/tSNE/Self_Correlation_Coefficients/
    mkdir $Main/Files/tSNE/Unitary_Vector_X/
    mkdir $Main/Files/tSNE/Unitary_Vector_Y/
    mkdir $Main/Files/tSNE/Transition_Probability/
    mkdir $Main/Files/tSNE/Transition_Probability/Sample_Level/
    mkdir $Main/Files/tSNE/Transition_Probability/Group_Level/
    tSNE="True"
else
    tSNE="False"
fi

if [ "$UMAPYesNo" = "yes" ]
then
    mkdir $Main/Fate_Maps/UMAP_2D/
    mkdir $Main/Files/UMAP_2D/
    mkdir $Main/Files/UMAP_2D/Delta_Embedding/
    mkdir $Main/Files/UMAP_2D/Self_Correlation_Coefficients/
    mkdir $Main/Files/UMAP_2D/Unitary_Vector_X/
    mkdir $Main/Files/UMAP_2D/Unitary_Vector_Y/
    mkdir $Main/Files/UMAP_2D/Transition_Probability/
    mkdir $Main/Files/UMAP_2D/Transition_Probability/Sample_Level/ 
    mkdir $Main/Files/UMAP_2D/Transition_Probability/Group_Level/
    UMAP="True"
else
    UMAP="False"
fi

if [ "$UMAP3DYesNo" = "yes" ]
then
    mkdir $Main/Fate_Maps/UMAP_3D/
    mkdir $Main/Files/UMAP_3D/
    mkdir $Main/Files/UMAP_3D/Delta_Embedding/
    mkdir $Main/Files/UMAP_3D/Self_Correlation_Coefficients/
    mkdir $Main/Files/UMAP_3D/Unitary_Vector_X/
    mkdir $Main/Files/UMAP_3D/Unitary_Vector_Y/
    mkdir $Main/Files/UMAP_3D/Unitary_Vector_Z/
    mkdir $Main/Files/UMAP_3D/Transition_Probability/
    mkdir $Main/Files/UMAP_3D/Transition_Probability/Sample_Level/
    mkdir $Main/Files/UMAP_3D/Transition_Probability/Group_Level/
    UMAP3D="True"
else
    UMAP3D="False"
fi

echo "----------------------------------------------------------"

cd "${MainExpressionFileLoc}/"

if [ "$RunMode" = "M1" ]
then
    metafiles=`ls Metadata_*`
    for eachfile in $metafiles;
    do
       #echo $eachfile
       sampname=${eachfile/"Metadata_"/}
       sampname=${sampname/".csv"/}
       splicedName="Spliced_"$sampname".csv"
       unsplicedName="Unspliced_"$sampname".csv"
       mkdir $Main/Fate_Maps/PCA/$sampname/
       if [ "$tSNEYesNo" = "yes" ]
       then
           mkdir $Main/Fate_Maps/tSNE/$sampname/
       fi
       if [ "$UMAPYesNo" = "yes" ]
       then
           mkdir $Main/Fate_Maps/UMAP_2D/$sampname/
       fi
       if [ "$UMAP3DYesNo" = "yes" ]
       then
           mkdir $Main/Fate_Maps/UMAP_3D/$sampname/
       fi
       mkdir $Main/Files/Velocity_Values/$sampname/
       mkdir $Main/Files/Predicted_Unspiced/$sampname/
       mkdir $Main/Files/Future_Spliced/$sampname/
       if [ "$UMAPYesNo" = "yes" ]
       then
           mkdir $Main/Files/UMAP_2D/Delta_Embedding/$sampname/
           mkdir $Main/Files/UMAP_2D/Self_Correlation_Coefficients/$sampname/
           mkdir $Main/Files/UMAP_2D/Unitary_Vector_X/$sampname/
           mkdir $Main/Files/UMAP_2D/Unitary_Vector_Y/$sampname/
           mkdir $Main/Files/UMAP_2D/Transition_Probability/Sample_Level/$sampname/
           mkdir $Main/Files/UMAP_2D/Transition_Probability/Group_Level/$sampname/
           if [ "$ConfidenceInterval" = "yes" ]
           then
               mkdir $Main/Files/UMAP_2D/Transition_Probability/Group_Level/$sampname/Confidence_Interval/
           fi
           if [ "$RandomSamp" = "yes" ]
           then
               mkdir $Main/Files/UMAP_2D/Transition_Probability/Group_Level/$sampname/Bootstrapping/
           fi
       fi
       if [ "$UMAP3DYesNo" = "yes" ]
       then
           mkdir $Main/Files/UMAP_3D/Delta_Embedding/$sampname/
           mkdir $Main/Files/UMAP_3D/Self_Correlation_Coefficients/$sampname/
           mkdir $Main/Files/UMAP_3D/Unitary_Vector_X/$sampname/
           mkdir $Main/Files/UMAP_3D/Unitary_Vector_Y/$sampname/
           mkdir $Main/Files/UMAP_3D/Unitary_Vector_Z/$sampname/
           mkdir $Main/Files/UMAP_3D/Transition_Probability/Sample_Level/$sampname/
           mkdir $Main/Files/UMAP_3D/Transition_Probability/Group_Level/$sampname/
           if [ "$ConfidenceInterval" = "yes" ]
           then
               mkdir $Main/Files/UMAP_3D/Transition_Probability/Group_Level/$sampname/Confidence_Interval/
           fi
           if [ "$RandomSamp" = "yes" ]
           then
               mkdir $Main/Files/UMAP_3D/Transition_Probability/Group_Level/$sampname/Bootstrapping/
           fi
       fi
       if [ "$tSNEYesNo" = "yes" ]
       then
           mkdir $Main/Files/tSNE/Delta_Embedding/$sampname/
           mkdir $Main/Files/tSNE/Self_Correlation_Coefficients/$sampname/
           mkdir $Main/Files/tSNE/Unitary_Vector_X/$sampname/
           mkdir $Main/Files/tSNE/Unitary_Vector_Y/$sampname/
           mkdir $Main/Files/tSNE/Transition_Probability/Sample_Level/$sampname/
           mkdir $Main/Files/tSNE/Transition_Probability/Group_Level/$sampname/
           if [ "$ConfidenceInterval" = "yes" ]
           then
               mkdir $Main/Files/tSNE/Transition_Probability/Group_Level/$sampname/Confidence_Interval/
           fi
           if [ "$RandomSamp" = "yes" ]
           then
               mkdir $Main/Files/tSNE/Transition_Probability/Group_Level/$sampname/Bootstrapping/
           fi
       fi
       if [ "$PCAYesNo" = "yes" ]
       then
           mkdir $Main/Files/PCA/Delta_Embedding/$sampname/
           mkdir $Main/Files/PCA/Self_Correlation_Coefficients/$sampname/
           mkdir $Main/Files/PCA/Unitary_Vector_X/$sampname/
           mkdir $Main/Files/PCA/Unitary_Vector_Y/$sampname/
           mkdir $Main/Files/PCA/Transition_Probability/Sample_Level/$sampname/
           mkdir $Main/Files/PCA/Transition_Probability/Group_Level/$sampname/
           if [ "$ThreeDPCA" = "yes" ]
           then
               mkdir $Main/Files/PCA/Unitary_Vector_Z/$sampname/
           fi
               if [ "$ConfidenceInterval" = "yes" ]
           then
               mkdir $Main/Files/PCA/Transition_Probability/Group_Level/$sampname/Confidence_Interval/
           fi
           if [ "$RandomSamp" = "yes" ]
           then
               mkdir $Main/Files/PCA/Transition_Probability/Group_Level/$sampname/Bootstrapping/
           fi
       fi
       cd "${MainExpressionFileLoc}/"   

       if [ $SecondMeta = "True" ]
       then
           Meta2Name="MetadataAlt_"$sampname".csv"
       else
           Meta2Name="Metadata_"$sampname".csv"
       fi
       if [ $RunMultiEmbed = "yes" ]
       then
           for i in $(seq $Embed1 $Embed2); do python3 -Wignore $CodeLocation/RNA_Velocity_Main.py $i "$wd" "$splicedName" "$unsplicedName" "${MainExpressionFileLoc}/${eachfile}" "${MainExpressionFileLoc}/$Meta2Name" "$PCAFS" "Results/" "$CCFSOption" "False" "$UMAPDistNum" "$UMAPDist" "$CodeLocation/" "$Typeof1" "$GIFRUN" "$SecondMeta" "$Typeof2" "$RS" "$CI" "$CIIteration" "$tSNE" "$UMAP" "$UMAP3D" "$PCANum" "$TPMin" "$TPMax" "$TPStep" "$sampname" "$ShapePoint" "$PCAN" "$ColorMap" "$PCAYesNo" "$ThreeDPCA" "$PredSecond" "$ProbabilityThreshold" "$AdjustPCANDS"; done
       else
           python3 -Wignore $CodeLocation/RNA_Velocity_Main.py $Embed1 "$wd" "$splicedName" "$unsplicedName" "${MainExpressionFileLoc}/${eachfile}" "${MainExpressionFileLoc}/$Meta2Name" "$PCAFS" "Results/" "$CCFSOption" "False" "$UMAPDistNum" "$UMAPDist" "$CodeLocation/" "$Typeof1" "$GIFRUN" "$SecondMeta" "$Typeof2" "$RS" "$CI" "$CIIteration" "$tSNE" "$UMAP" "$UMAP3D" "$PCANum" "$TPMin" "$TPMax" "$TPStep" "$sampname" "$ShapePoint" "$PCAN" "$ColorMap" "$PCAYesNo" "$ThreeDPCA" "$PredSecond" "$ProbabilityThreshold" "$AdjustPCANDS"
       fi
    done
else
    sampname=$CommonSubstring
    metafile="Metadata_"$sampname".csv"
    splicedName="Spliced_"$sampname".csv"
    unsplicedName="Unspliced_"$sampname".csv"
    mkdir $Main/Fate_Maps/PCA/$sampname/
    if [ "$tSNEYesNo" = "yes" ]
    then
        mkdir $Main/Fate_Maps/tSNE/$sampname/
    fi
    if [ "$UMAPYesNo" = "yes" ]
    then
        mkdir $Main/Fate_Maps/UMAP_2D/$sampname/
    fi
    if [ "$UMAP3DYesNo" = "yes" ]
    then
        mkdir $Main/Fate_Maps/UMAP_3D/$sampname/
    fi
    mkdir $Main/Files/Velocity_Values/$sampname/
    mkdir $Main/Files/Predicted_Unspiced/$sampname/
    mkdir $Main/Files/Future_Spliced/$sampname/
    if [ "$UMAPYesNo" = "yes" ]
    then
        mkdir $Main/Files/UMAP_2D/Delta_Embedding/$sampname/
        mkdir $Main/Files/UMAP_2D/Self_Correlation_Coefficients/$sampname/
        mkdir $Main/Files/UMAP_2D/Unitary_Vector_X/$sampname/
        mkdir $Main/Files/UMAP_2D/Unitary_Vector_Y/$sampname/
        mkdir $Main/Files/UMAP_2D/Transition_Probability/Sample_Level/$sampname/
        mkdir $Main/Files/UMAP_2D/Transition_Probability/Group_Level/$sampname/
        if [ "$ConfidenceInterval" = "yes" ]
        then
            mkdir $Main/Files/UMAP_2D/Transition_Probability/Group_Level/$sampname/Confidence_Interval/
        fi
        if [ "$RandomSamp" = "yes" ]
        then
            mkdir $Main/Files/UMAP_2D/Transition_Probability/Group_Level/$sampname/Bootstrapping/
        fi
    fi
    if [ "$UMAP3DYesNo" = "yes" ]
    then
        mkdir $Main/Files/UMAP_3D/Delta_Embedding/$sampname/
        mkdir $Main/Files/UMAP_3D/Self_Correlation_Coefficients/$sampname/
        mkdir $Main/Files/UMAP_3D/Unitary_Vector_X/$sampname/
        mkdir $Main/Files/UMAP_3D/Unitary_Vector_Y/$sampname/
        mkdir $Main/Files/UMAP_3D/Unitary_Vector_Z/$sampname/
        mkdir $Main/Files/UMAP_3D/Transition_Probability/Sample_Level/$sampname/
        mkdir $Main/Files/UMAP_3D/Transition_Probability/Group_Level/$sampname/
        if [ "$ConfidenceInterval" = "yes" ]
        then
            mkdir $Main/Files/UMAP_3D/Transition_Probability/Group_Level/$sampname/Confidence_Interval/
        fi
        if [ "$RandomSamp" = "yes" ]
        then
            mkdir $Main/Files/UMAP_3D/Transition_Probability/Group_Level/$sampname/Bootstrapping/
        fi
    fi
    if [ "$tSNEYesNo" = "yes" ]
    then
        mkdir $Main/Files/tSNE/Delta_Embedding/$sampname/
        mkdir $Main/Files/tSNE/Self_Correlation_Coefficients/$sampname/
        mkdir $Main/Files/tSNE/Unitary_Vector_X/$sampname/
        mkdir $Main/Files/tSNE/Unitary_Vector_Y/$sampname/
        mkdir $Main/Files/tSNE/Transition_Probability/Sample_Level/$sampname/
        mkdir $Main/Files/tSNE/Transition_Probability/Group_Level/$sampname/
        if [ "$ConfidenceInterval" = "yes" ]
        then
            mkdir $Main/Files/tSNE/Transition_Probability/Group_Level/$sampname/Confidence_Interval/
        fi
        if [ "$RandomSamp" = "yes" ]
        then
            mkdir $Main/Files/tSNE/Transition_Probability/Group_Level/$sampname/Bootstrapping/
        fi
    fi
    if [ "$PCAYesNo" = "yes" ]
    then
        mkdir $Main/Files/PCA/Delta_Embedding/$sampname/
        mkdir $Main/Files/PCA/Self_Correlation_Coefficients/$sampname/
        mkdir $Main/Files/PCA/Unitary_Vector_X/$sampname/
        mkdir $Main/Files/PCA/Unitary_Vector_Y/$sampname/
        mkdir $Main/Files/PCA/Transition_Probability/Sample_Level/$sampname/
        mkdir $Main/Files/PCA/Transition_Probability/Group_Level/$sampname/
        if [ "$ThreeDPCA" = "yes" ]
        then
            mkdir $Main/Files/PCA/Unitary_Vector_Z/$sampname/
        fi
        if [ "$ConfidenceInterval" = "yes" ]
        then
            mkdir $Main/Files/PCA/Transition_Probability/Group_Level/$sampname/Confidence_Interval/
        fi
        if [ "$RandomSamp" = "yes" ]
        then
            mkdir $Main/Files/PCA/Transition_Probability/Group_Level/$sampname/Bootstrapping/
        fi
    fi
    cd "${MainExpressionFileLoc}/"
    if [ "$SecondMeta" = "True" ]
    then
        Meta2Name="MetadataAlt_"$sampname".csv"
    else
        Meta2Name="Metadata_"$sampname".csv"
    fi
    if [ "$RunMultiEmbed" = "yes" ]
    then
        for i in $(seq $Embed1 $Embed2); do python3 -Wignore $CodeLocation/RNA_Velocity_Main.py $i "$wd" "$splicedName" "$unsplicedName" "${MainExpressionFileLoc}/$metafile" "${MainExpressionFileLoc}/$Meta2Name" "$PCAFS" "Results/" "$CCFSOption" "False" "$UMAPDistNum" "$UMAPDist" "$CodeLocation/" "$Typeof1" "$GIFRUN" "$SecondMeta" "$Typeof2" "$RS" "$CI" "$CIIteration" "$tSNE" "$UMAP" "$UMAP3D" "$PCANum" "$TPMin" "$TPMax" "$TPStep" "$sampname" "$ShapePoint" "$PCAN" "$ColorMap" "$PCAYesNo" "$ThreeDPCA" "$PredSecond" "$ProbabilityThreshold" "$AdjustPCANDS"; done
    else
        python3 -Wignore $CodeLocation/RNA_Velocity_Main.py $Embed1 "$wd" "$splicedName" "$unsplicedName" "${MainExpressionFileLoc}/$metafile" "${MainExpressionFileLoc}/$Meta2Name" "$PCAFS" "Results/" "$CCFSOption" "False" "$UMAPDistNum" "$UMAPDist" "$CodeLocation/" "$Typeof1" "$GIFRUN" "$SecondMeta" "$Typeof2" "$RS" "$CI" "$CIIteration" "$tSNE" "$UMAP" "$UMAP3D" "$PCANum" "$TPMin" "$TPMax" "$TPStep" "$sampname" "$ShapePoint" "$PCAN" "$ColorMap" "$PCAYesNo" "$ThreeDPCA" "$PredSecond" "$ProbabilityThreshold" "$AdjustPCANDS"
    fi

fi

if [ "$SensSpec" = "yes" ] || [ "$SensSpec2" = "yes" ]
then
    mkdir $Main/Summary/
    echo "Calculating sensitivity and specificity values..."
fi

if [ "$SensSpec" = "yes" ]
then
    output="$Main/Summary/"
    if [ "$tSNE" = "True" ]
    then
        resfilers="$Main/Files/tSNE/Transition_Probability/Group_Level/$sampname/"
        python3 $CodeLocation/SensSpecScript.py "$resfilers" "$sampname" "${MainExpressionFileLoc}/Metadata_$sampname.csv" "tSNE" "$output" "$Feature1" "$Feature2" "tSNE" "SampleGroup1" "$CodeLocation/" "$ProbabilityThreshold"
    fi
    if [ "$UMAP" = "True" ]
    then
        resfilers2="$Main/Files/UMAP_2D/Transition_Probability/Group_Level/$sampname/"
        python3 $CodeLocation/SensSpecScript.py "$resfilers2" "$sampname" "${MainExpressionFileLoc}/Metadata_$sampname.csv" "UMAP_2D" "$output" "$Feature1" "$Feature2" "UMAP" "SampleGroup1" "$CodeLocation/" "$ProbabilityThreshold"
    fi
    if [ "$UMAP3D" = "True" ]
    then
        resfilers3="$Main/Files/UMAP_3D/Transition_Probability/Group_Level/$sampname/"
        python3 $CodeLocation/SensSpecScript.py "$resfilers3" "$sampname" "${MainExpressionFileLoc}/Metadata_$sampname.csv" "UMAP_3D" "$output" "$Feature1" "$Feature2" "UMAP3D" "SampleGroup1" "$CodeLocation/" "$ProbabilityThreshold"
    fi
    if [ "$PCAYesNo" = "yes" ]
    then
        resfilers4="$Main/Files/PCA/Transition_Probability/Group_Level/$sampname/"
        python3 $CodeLocation/SensSpecScript.py "$resfilers4" "$sampname" "${MainExpressionFileLoc}/Metadata_$sampname.csv" "PCA" "$output" "$Feature1" "$Feature2" "PCA_2D" "SampleGroup1" "$CodeLocation/" "$ProbabilityThreshold"
        if [ "$ThreeDPCA" = "yes" ]
        then
            python3 $CodeLocation/SensSpecScript.py "$resfilers4" "$sampname" "${MainExpressionFileLoc}/Metadata_$sampname.csv" "PCA" "$output" "$Feature1" "$Feature2" "PCA_3D" "SampleGroup1" "$CodeLocation/" "$ProbabilityThreshold"
        fi
    fi
fi


if [ "$SensSpec2" = "yes" ]
then
    output="$Main/Summary/"
    if [ "$tSNE" = "True" ]
    then
        resfilers5="$Main/Files/tSNE/Transition_Probability/Group_Level/$sampname/"
        python3 $CodeLocation/SensSpecScript.py "$resfilers5" "$sampname" "${MainExpressionFileLoc}/MetadataAlt_$sampname.csv" "tSNE" "$output" "$Feature21" "$Feature22" "tSNE" "SampleGroup2" "$CodeLocation/" "$ProbabilityThreshold"
    fi
    if [ "$UMAP" = "True" ]
    then
        resfilers6="$Main/Files/UMAP_2D/Transition_Probability/Group_Level/$sampname/"
        python3 $CodeLocation/SensSpecScript.py "$resfilers6" "$sampname" "${MainExpressionFileLoc}/MetadataAlt_$sampname.csv" "UMAP_2D" "$output" "$Feature21" "$Feature22" "UMAP" "SampleGroup2" "$CodeLocation/" "$ProbabilityThreshold"
    fi
    if [ "$UMAP3D" = "True" ]
    then
        resfilers7="$Main/Files/UMAP_3D/Transition_Probability/Group_Level/$sampname/"
        python3 $CodeLocation/SensSpecScript.py "$resfilers7" "$sampname" "${MainExpressionFileLoc}/MetadataAlt_$sampname.csv" "UMAP_3D" "$output" "$Feature21" "$Feature22" "UMAP3D" "SampleGroup2" "$CodeLocation/" "$ProbabilityThreshold"
    fi
    if [ "$PCAYesNo" = "yes" ]
    then
        resfilers8="$Main/Files/PCA/Transition_Probability/Group_Level/$sampname/"
        python3 $CodeLocation/SensSpecScript.py "$resfilers8" "$sampname" "${MainExpressionFileLoc}/MetadataAlt_$sampname.csv" "PCA" "$output" "$Feature21" "$Feature22" "PCA_2D" "SampleGroup2" "$CodeLocation/" "$ProbabilityThreshold"
        if [ "$ThreeDPCA" = "yes" ]
        then
            python3 $CodeLocation/SensSpecScript.py "$resfilers8" "$sampname" "${MainExpressionFileLoc}/MetadataAlt_$sampname.csv" "PCA" "$output" "$Feature21" "$Feature22" "PCA_3D" "SampleGroup2" "$CodeLocation/" "$ProbabilityThreshold"
        fi
    fi
fi


echo "Congratulations - the analysis is finished!"

echo "----------------------------------------------------------"

#parameters dor RNA velocity python code:
#1 is the number of neighbours/perplexity
#2 loction of everything
#3 spliced file
#4 unspliced file
#5 main metadata file with samples as first columan and Group as second ("Group" must be second column name) - used for colouts
#6 second metadata file of additional info - used for markers and/or as additional thing to look at trabsition probabilities over
#7 PCA variation based feature selection: NoPCAVar (do not run) or PCAVar (run this analysis)
#8 location of output dir (to be added to the end of the second parameter
#9 Do basic correlation coefficient based feature selection? NoCoorelationCoeff (do not run) or CoorelationCoeff (run this)
#10 apply scaling? (for vlm.calculate_embedding_shift* functions) False (do not do this) or True (do this)
#11 distance parameter for UMAP (currently set to the default)
#12 metric method used for UMAP parameter (currently set the default)
#13 location of the additional (R) code required for the transition probability work
#14 the "type" of (grouptype) first metadata i.e. are these "strings" (classes) or "integers" (numeric timepoints)
#15 generate gif versions of UMAP3D, transiton probability scatterplots and 3D PCA? True (run this) or False (do not run this)
#16 true is "manyfeatures" - basically if your second ematadat file has useful classes to predict on then put as "True"
#17 "strings": type of second feature is the final paramater
#18 true is Randomsampling
#19 confidence interval calculation? (using first metadata: 7)) - "True" or "False"
#20 number of iterations for CI calculation - always an integer
#21 run tsne analsysi? "True" or "False"
#22 run UMAP analysis? "True" or "False" - second command line argument for the shell script
#23 run 3d umap? "True" or "False" - if false then it will just run the 2d umap, has to be used in conjuncture witha true argument 24 - third commend line argument for the shell script
#24 number of principal components to use for the pca! def: "3"  must be an integer defines how many PCs are given to tSNE
#25 beginning of the range for the transtion probability cloest number of neighbours value
#26 the end of a range for the transtion probability cloest number of neighbours value i.e. "2" at #27 and "11" at #28 would generate values between 2-10
#27 the skip i.e. betwen 27 and 28 - 2 would be 2,4,6, 1 would be 2,3,4
#28 specific output folder name
#29 Whether to shape the points by the metadata
#30 number of pca components to plot i.e. first "2" or first "3"
#31 colour map
#32 fate maps pca plots
#33 3d pca fate maps
#34 use the second metadata feature for prediction - yes or no
#35 user threshold for prediction


