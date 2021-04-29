// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Engine/Texture2D.h"
#include "EditorUtilityWidget.h"
#include "LanGenEditorUtilityWidget.generated.h"

/**
 * 
 */
UCLASS(BlueprintType)
class LANDSCAPEGENERATOR_API ULanGenEditorUtilityWidget : public UEditorUtilityWidget
{
	GENERATED_BODY()
	
public:
	ULanGenEditorUtilityWidget();
	~ULanGenEditorUtilityWidget();

	UFUNCTION(BlueprintCallable)
		void ConLog(FString text);
	UFUNCTION(BlueprintCallable)
		void ScrLog(FString text);
	UFUNCTION(BlueprintCallable)
		UTexture2D* GenerateTexture(int x, int y);

	UFUNCTION(BlueprintImplementableEvent)
		void GenerateFunction(const int& xCoord, const int& yCoord);

	UFUNCTION(BlueprintCallable)
		float NoiseTest(int xCoord, int yCoord);
};
