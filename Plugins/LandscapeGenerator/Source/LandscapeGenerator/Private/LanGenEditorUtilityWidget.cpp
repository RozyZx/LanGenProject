// Fill out your copyright notice in the Description page of Project Settings.


#include "LanGenEditorUtilityWidget.h"
#include "Engine/Texture2D.h"
#include <Runtime/Landscape/Classes/Landscape.h>
#include <Landscape.h>

#define SCR_LOG(x, ...) if(GEngine){GEngine->AddOnScreenDebugMessage(-1, 15.0f, FColor::Red, FString::Printf(TEXT(x), __VA_ARGS__));}
#define CON_LOG(x, ...) UE_LOG(LogTemp, Warning, TEXT(x), __VA_ARGS__);

ULanGenEditorUtilityWidget::ULanGenEditorUtilityWidget() {}
ULanGenEditorUtilityWidget::~ULanGenEditorUtilityWidget() {}

void ULanGenEditorUtilityWidget::ConLog(FString text) { CON_LOG("%s", *text); }

void ULanGenEditorUtilityWidget::ScrLog(FString text) { SCR_LOG("%s", *text); }

UTexture2D* ULanGenEditorUtilityWidget::GenerateTexture(int x, int y)
{
	UTexture2D* texture = UTexture2D::CreateTransient(x, y, PF_B8G8R8A8);
	TArray<int8> pixelData;						//generate red only
	pixelData.Init(0, x * y * 4);

	/*for (int32 xx = 0; xx < x; ++xx) {
		for (int32 yy = 0; yy < y; ++yy) {
			//GenerationFunctions(xx, yy)*255
			pixelData.Insert(GenerateFunction(xx, yy) * 255, 0 + (y * yy) + (x * xx));
		}
	}*/

	void* TextureData = texture->PlatformData->Mips[0].BulkData.Lock(LOCK_READ_WRITE);
	void* data = &pixelData;
	FMemory::Memcpy(TextureData, data, 4 * x * y);
	texture->PlatformData->Mips[0].BulkData.Unlock();
	texture->UpdateResource();

	return texture;
}

float ULanGenEditorUtilityWidget::NoiseTest(int xCoord, int yCoord)
{
	return 0.0;
}