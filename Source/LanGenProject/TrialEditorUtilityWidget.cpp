// Fill out your copyright notice in the Description page of Project Settings.


#include "TrialEditorUtilityWidget.h"
#include "Engine/Texture2D.h"
#include <Runtime/Landscape/Classes/Landscape.h>
#include <Landscape.h>
#include "Math/UnrealMathUtility.h"
#include "ImageUtils.h"
#include "Math/Color.h"

#define SCR_LOG(x, ...) if(GEngine){GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Red, FString::Printf(TEXT(x), __VA_ARGS__));}
#define CON_LOG(x, ...) UE_LOG(LogTemp, Warning, TEXT(x), __VA_ARGS__);

void UTrialEditorUtilityWidget::ConLog(FString text) { CON_LOG("%s", *text); }

void UTrialEditorUtilityWidget::ScrLog(FString text) { SCR_LOG("%s", *text); }

UTexture2D* UTrialEditorUtilityWidget::GenerateTexture(int x, int y, int tileX, int tileY, int generateParam)
{
	TArray<FColor> colorData;
	int index = 0;
	FVector2D location;
	GenerationSetup();

	switch (generateParam) {
	case 0: colorData = GenerateElevationMap(x, y); break;
	case 1: colorData.Init(FColor::Black, x * y); break; //noise only
	case 2: colorData = GenerateElevationMap(x, y); goto create; //elevation only
	}
	ConLog("noise creation");
	for (int i = 0; i < x; ++i) {
		for (int j = 0; j < y; ++j) {
			index = (i * y) + j;
			location.X = (float) i / tileX;
			location.Y = (float) j / tileY;
			colorData[index].R = GenerateFunction(location, colorData[index].R);
		}
	}

	create:
	ConLog(FString::FromInt(colorData.Num()));
	//colorData.Init(FColor::Black, x * y);
	FCreateTexture2DParameters textureParam;
	UTexture2D* texture = FImageUtils::CreateTexture2D(x, y, colorData, this, TEXT("heightTexture"), RF_NoFlags, textureParam);

	return texture;
}

uint8 UTrialEditorUtilityWidget::MapTo8Bit(float in, float min, float max) { return 255 * ((in - min) / (max - min)); }

int UTrialEditorUtilityWidget::MapFloatToInt(float in, float inMin, float inMax, int min, int max)
{
	return FMath::RoundToInt((max - min) * (in / (inMax - inMin)));
}

float UTrialEditorUtilityWidget::RandomFloat() { return FMath::RandRange((float)-1.0, (float)1.0); }

void UTrialEditorUtilityWidget::Debug(int x, int y, int tileX, int tileY)
{
	TArray<FColor> colorData = GenerateElevationMap(x, y);
	int index = 0;
	FVector2D location;
	GenerationSetup();

	ConLog(FString::FromInt(colorData.Num()));
}
