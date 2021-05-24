/*
* Bresenham converted from https://gist.github.com/polyamide/424a9fc8507436883a911e6f1d79d9d3#file-bresenham-js-L1
*
* Edited and converted to Unreal compatible c++ by Fachrurrozy Muhammad
*/

#include "LanGenElevation.h"
#include "Math/UnrealMathUtility.h"
#include "Misc/Char.h"
#include "Misc/DefaultValueHelper.h"

#define CON_LOG(x, ...) UE_LOG(LogTemp, Warning, TEXT(x), __VA_ARGS__);

void ULanGenElevation::ResetSeed()
{
	p.Empty();
	p = P_BASE;
	p.Append(P_BASE);
}

void ULanGenElevation::Init(int32 in, int x, int y)
{
    lanX = x;
    lanY = y;
	seed = in;
	p.Empty();
	p = P_BASE;
	randomEngine = FRandomStream(seed);
	TArray<int> tempArray = P_BASE;
	Shuffle(tempArray);
	p = tempArray;
	p.Append(tempArray);
}

FVector2D ULanGenElevation::RandomizeCoord(float percentageSafeZone)
{
    float half = (1 - percentageSafeZone) / 2;
    int x = randomEngine.RandRange(half * lanX, (percentageSafeZone + half) * lanX),
        y = randomEngine.RandRange(half * lanY, (percentageSafeZone + half) * lanY);
    return FVector2D(
        randomEngine.RandRange(half * lanX, (percentageSafeZone + half) * lanX),
        randomEngine.RandRange(half * lanY, (percentageSafeZone + half) * lanY)
    );
}

FVector2D ULanGenElevation::CenterTerrainCoord() { return FVector2D(lanX / 2, lanY / 2); }

TArray<FColor> ULanGenElevation::GenerateGraph(
    FVector2D startingPosition, FString rule, FString axiom,
    int ruleLoop, int lineLength, int minAngle,
    int maxAngle, int radius, int peak,
    float skew, int fillDegree, float topBlend,
    int disLoop, float disSmooth, int startHeight
)
{
    TArray<coord> branchRootStack;
    TArray<coord> currentLine;
    coord* currentCoord;
    FString grammar;
    bool isRandomAngle = minAngle != maxAngle;
    int peakIndex = 0;

    init = FColor(startHeight, 0, 0);
    texture.Init(init, lanX * lanY);
    detailTexture.Init(FColor::Black, lanX * lanY);

    // L-System
    RuleSetup(rule);
    grammar = RuleApply(axiom, ruleLoop);
    //CON_LOG("%s", *grammar);

    // set starting position
    currentLine.Add(coord(startingPosition.X, startingPosition.Y, 0));
    CON_LOG("Start: (%d, %d)", currentLine[0].x, currentLine[0].y);

    // create array of target coord
    for (TCHAR i : grammar) {
        currentCoord = &currentLine[currentLine.Num() - 1];
        switch (i) {
        case 'F': Bresenham(currentLine, lineLength); break;
        case 'P': peakIndex = currentLine.Num() - 1; break;
        case 'L': /*main lowest point*/
            if (currentLine.Num() > 2) {
                MidpointDisplacement(currentLine, peak, peakIndex, peak / 2, disLoop, disSmooth);
                GradientSingleMain(currentLine, peak, radius, skew, fillDegree, topBlend, false);
            }
            break;
        case 'D': Bresenham(currentLine, lineLength); break;
        case 'E': /*draw detail*/
            if (currentLine.Num() > 1) {
                for (coord& curCoord : currentLine) curCoord.height = peak * 0.1;
                GradientSingleMain(currentLine, currentCoord->height, 0.1 * radius, 0, 180, 0.5 * topBlend, false, true);
            }
            break;
        case '+': currentCoord->AddTheta(isRandomAngle ? randomEngine.RandRange(minAngle, maxAngle) : minAngle); break;
        case '-': currentCoord->AddTheta(-(isRandomAngle ? randomEngine.RandRange(minAngle, maxAngle) : minAngle)); break;
        case '[': branchRootStack.Add(*currentCoord); break;
        case ']': // run on line ends
            currentLine.Empty();
            currentLine.Add(
                branchRootStack[branchRootStack.Num() - 1]
            );
            branchRootStack.RemoveAt(branchRootStack.Num() - 1);
            break;
        }
    }

    for (int i = 0; i < texture.Num(); ++i) {
        texture[i].R += detailTexture[i].R;
    }

    return texture;

}

TArray<FColor> ULanGenElevation::CombineTexture(TArray<FColor> texture1, TArray<FColor> texture2)
{
    if (texture2.Num() != texture1.Num()) return texture1;
    for (int i = 0; i < texture1.Num(); ++i) texture1[i].R = texture1[i].R >= texture2[i].R ? texture1[i].R : texture2[i].R;
    return texture1;
}

// L-System rule setup
void ULanGenElevation::RuleSetup(FString in)
{
    /* in = F{[F]F:25,-F:25,+F:25,FF:25} */
    rules.Empty();
    TArray<FString> inRules, tooRule;
    rule inRule;
    FString from, too, to, probString;
    int prob;

    in.ParseIntoArray(inRules, TEXT("}"));
    for (FString i : inRules) {
        i.Split("{", &from, &too);
        inRule.from = from;
        too.ParseIntoArray(tooRule, TEXT(","));
        for (FString j : tooRule) {
            j.Split(":", &to, &probString);
            inRule.to.Add(to);
            FDefaultValueHelper::ParseInt(probString, prob);
            inRule.prob.Add(prob);
        }
        rules.Add(inRule);
        inRule.from.Empty();
        inRule.to.Empty();
        inRule.prob.Empty();
    }

    // debugging purpose
    /*for (rule i : rules) {
        CON_LOG("from: %s", *i.from);
        for (FString j : i.to) CON_LOG("to: %s", *j);
    }*/
}

// apply L-System rule
FString ULanGenElevation::RuleApply(FString axiom, int loop)
{
    FString currentLstring, last;
    int lCount;
    for (int i = 0; i < loop; ++i) {
        lCount = 0;
        last = "";
        currentLstring = "";
        //check for each char in axiom
        for (TCHAR j : axiom) {
            if (lCount < 2) {
                if (j == 'L') ++lCount;
                for (int k = 0; k < rules.Num(); ++k) {
                    if (j == rules[k].from[0]) {
                        currentLstring.Append(RandomizeRule(k, last));
                        break;
                    }
                }
                currentLstring.AppendChar(j);

                if (lCount == 2) { // reset variable after use
                    currentLstring.Append(last);
                    last = "";
                }
            }
            else {
                // rule check
                for (int k = 0; k < rules.Num(); ++k) {
                    if (j == rules[k].from[0]) {
                        currentLstring.Append(RandomizeRule(k, last));
                        break;
                    }
                }
                // if nothing match
                currentLstring.AppendChar(j);

                if (j == 'E') {
                    currentLstring.Append(last);
                    last = "";
                }
            }
        }
        axiom = currentLstring;
        if (axiom.Len() > 1000000000) break; // prevent editor from crashing due string length limit
    }
    return axiom;
}

void ULanGenElevation::Shuffle(TArray<int>& inArr)
{
    for (int i = 0; i < inArr.Num(); ++i) inArr.Swap(i, randomEngine.RandRange(0, p.Num() - 1));
}

// decide which symbol to replace to based on rule
FString ULanGenElevation::RandomizeRule(int ruleIndex, FString& last)
{
    int randomNumber = randomEngine.RandRange(1, 100);
    int currentRequirement = 0;
    FString res;
    bool addLast = false;
    for (int i = 0; i < rules[ruleIndex].to.Num(); ++i) {
        currentRequirement += rules[ruleIndex].prob[i];
        if (randomNumber <= currentRequirement) {
            for (TCHAR j : rules[ruleIndex].to[i]) {
                if (j == ']') addLast = true;

                if (addLast) last.AppendChar(j);
                else res.AppendChar(j);
            }
            return res;
        }
    }
    res = rules[ruleIndex].from;
    return res;
}

/*
FString ULanGenElevation::RandomizeRule(int ruleIndex, FString& last)
{
    int randomNumber = randomEngine.RandRange(1, 100);
    int currentRequirement = 0;
    for (int i = 0; i < rules[ruleIndex].to.Num(); ++i) {
        currentRequirement += rules[ruleIndex].prob[i];
        if (randomNumber <= currentRequirement) {
            return rules[ruleIndex].to[i];
        }
    }
    return rules[ruleIndex].from;
}
*/

// get all 'hit' coord
void ULanGenElevation::Bresenham(TArray<coord>& currentLine, int lineLength)
{
    TArray<coord> res;
    coord currentCoord = currentLine[currentLine.Num() - 1];

    int x0 = currentCoord.x;
    int y0 = currentCoord.y;
    int x = x0 + FMath::Sin(DegreeToRad(currentCoord.theta)) * lineLength;
    int y = y0 + FMath::Cos(DegreeToRad(currentCoord.theta)) * lineLength;

    int dx = FMath::Abs(x - x0);
    int dy = FMath::Abs(y - y0);
    int sx = (x0 < x) ? 1 : -1;
    int sy = (y0 < y) ? 1 : -1;
    int err = dx - dy;
    int e2;

    while (!((x0 == x) && (y0 == y))) {
        e2 = err << 1;
        if (e2 > -dy) {
            err -= dy;
            x0 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y0 += sy;
        }
        currentLine.Add(coord(x0, y0, currentCoord.theta));
        //index = x0 * lanY + y0;
    }
}

// generate height for the ridge line
void ULanGenElevation::MidpointDisplacement(TArray<coord>& currentLineCoord, int peak, int peakIndex, int displacement, int loop, float smooth)
{
    TArray<midPoint> oldIndexes, newIndexes;
    midPoint mid;
    int currentDisplacement = displacement,
        randomModifier = 0;
    float modifier = FMath::Pow(2, -smooth), linearM;

    // first half
    if (FMath::Pow(2, loop) > peakIndex) loop = FMath::Log2(peakIndex + 1);

    oldIndexes.Add(midPoint(0, 0));
    oldIndexes.Add(midPoint(peakIndex, peak));

    //if (FMath::Pow(2, loop) >= currentLineCoord.Num()) {
    //    loop = FMath::Log2(currentLineCoord.Num());
    //}
    //if (currentLineCoord[0].isInRange(lanX, lanY)) {
    //    oldIndexes.Add(midPoint(0,
    //        texture[currentLineCoord[0].index(lanY)].R - init.R));
    //}
    //else oldIndexes.Add(midPoint(0, 0));
    //if (currentLineCoord[currentLineCoord.Num() - 1].isInRange(lanX, lanY)) {
    //    oldIndexes.Add(midPoint(currentLineCoord.Num() - 1,
    //        texture[currentLineCoord[currentLineCoord.Num() - 1].index(lanY)].R - init.R));
    //}
    //else oldIndexes.Add(midPoint(currentLineCoord.Num() - 1, 0));
    //// if branch
    //if (oldIndexes[0].height > 0) currentDisplacement *= modifier;

    // fill height between coord
    linearM = LinearM(peak, peakIndex + 1, 1);
    for (int j = 0; j < oldIndexes.Num() - 1; ++j) 
        for (int k = oldIndexes[j].index; k < oldIndexes[j + 1].index; ++k) currentLineCoord[k].height = Linear(linearM, k - peakIndex, peak);
    // midpoint
    for (int i = 0; i < loop; ++i) {
        for (int j = 0; j < oldIndexes.Num() - 1; ++j) {
            mid.index = (int)((oldIndexes[j].index + oldIndexes[j + 1].index) / 2);
            newIndexes.Add(oldIndexes[j]);
            if (mid.index != oldIndexes[j].index) {
                currentDisplacement *= modifier;
                /*if (currentLineCoord[mid.index].height < currentDisplacement) randomModifier = 1;
                else randomModifier = (randomEngine.RandRange(0, 1) == 0 ? 1 : -1);*/
                mid.height = currentLineCoord[mid.index].height + currentDisplacement * (randomEngine.RandRange(0,1) == 0 ? -1 : 1);
                newIndexes.Add(mid);
            }
        }
        newIndexes.Add(oldIndexes[oldIndexes.Num() - 1]);
        oldIndexes = newIndexes;
        newIndexes.Empty();
        // fill height between coord
        for (int j = 0; j < oldIndexes.Num() - 1; ++j) {
            for (int k = oldIndexes[j].index; k < oldIndexes[j + 1].index; ++k) {
                currentLineCoord[k].height =
                    (
                        ((float)(k - oldIndexes[j].index) / (oldIndexes[j + 1].index - oldIndexes[j].index)) *
                        (oldIndexes[j + 1].height - oldIndexes[j].height)) +
                    oldIndexes[j].height;
            }
        }
    }
    
    // second half
    oldIndexes.Empty();
    newIndexes.Empty();
    currentDisplacement = displacement;
    if (FMath::Pow(2, loop) > currentLineCoord.Num() - peakIndex) loop = FMath::Log2(peakIndex + 1);
    oldIndexes.Add(midPoint(peakIndex, peak));
    oldIndexes.Add(midPoint(currentLineCoord.Num() - 1, 0));
    // fill height between coord
    linearM = LinearM(peak, currentLineCoord.Num() - 1 - peakIndex);
    for (int j = 0; j < oldIndexes.Num() - 1; ++j) 
        for (int k = oldIndexes[j].index; k < oldIndexes[j + 1].index; ++k) currentLineCoord[k].height = Linear(linearM, k - peakIndex, peak);
    // midpoint
    for (int i = 0; i < loop; ++i) {
        for (int j = 0; j < oldIndexes.Num() - 1; ++j) {
            mid.index = (int)((oldIndexes[j].index + oldIndexes[j + 1].index) / 2);
            newIndexes.Add(oldIndexes[j]);
            if (mid.index != oldIndexes[j].index) {
                currentDisplacement *= modifier;
                mid.height = currentLineCoord[mid.index].height + currentDisplacement * (randomEngine.RandRange(0, 1) == 0 ? -1 : 1);
                newIndexes.Add(mid);
            }
        }
        newIndexes.Add(oldIndexes[oldIndexes.Num() - 1]);
        oldIndexes = newIndexes;
        newIndexes.Empty();
        // fill height between coord
        for (int j = 0; j < oldIndexes.Num() - 1; ++j) {
            for (int k = oldIndexes[j].index; k < oldIndexes[j + 1].index; ++k) {
                currentLineCoord[k].height =
                    (
                        ((float)(k - oldIndexes[j].index) / (oldIndexes[j + 1].index - oldIndexes[j].index)) *
                        (oldIndexes[j + 1].height - oldIndexes[j].height)) +
                    oldIndexes[j].height;
            }
        }
    }
}

void ULanGenElevation::GradientSingleMain(TArray<coord>& curLine, int peak, int radius, float skew, int fillDegree, float topBlend, bool calcHeight, bool isAdd)
{
    // main line height
    int mainLineRadius = curLine.Num() / 2,
        mainLineParaRadius = mainLineRadius * 0.75,
        /*mainLineExpRadius = mainLineRadius * 0.25,*/
        blendHeight, blendOffset;
    float paraA = ParabolaA(peak, mainLineParaRadius);
    float expA = ExponentDecayA(peak, mainLineParaRadius); // use para radius for smoother result
    blendOffset = ParabolaX(paraA, peak, peak / 4);
    blendHeight = (peak / 4) / expA;

    // ridgeline height
    if (calcHeight) {
        for (int i = 0; i < curLine.Num(); ++i) {
            if (i - mainLineRadius < -blendOffset) {
                // left blend
                curLine[i].height = ExponentDecay(expA, blendHeight, i, mainLineRadius - blendOffset, -1);;
            }
            else if (i - mainLineRadius > blendOffset) {
                // right blend
                curLine[i].height = ExponentDecay(expA, blendHeight, i, mainLineRadius + blendOffset);
            }
            else {
                curLine[i].height = Parabola(paraA, peak, i, mainLineRadius);
            }
            // crash *memory usage* prevention
            if (curLine[i].height < 0) curLine[i].height = 0;
            //curLine[i].height = peak;
        }
    }
    // gradient; peak -> x; radius -> y; y = mx; m = y / x
    int eachRadius;
    float m = (float) radius / peak;

    for (coord curCoord : curLine) {
        eachRadius = m * curCoord.height;
        //CON_LOG("m calc: %f \teach rad: %d, \theight: %d", (float) radius / peak, eachRadius, curCoord.height);
        GradientSingleMainHelper(curCoord, eachRadius, skew, fillDegree, topBlend, isAdd);
    }
}

void ULanGenElevation::GradientSingleMainHelper(coord curCoord, int radius, float skew, int fillDegree, float topBlend, bool isAdd)
{
    int leftRadius = radius * ((-skew) + 1),
        rightRadius = radius * (skew + 1),
        curIndex = 0, tempTheta = 0,
        blendLeftHeight = 0, blendRightHeight = 0,
        blendLeftOffset = 0, blendRightOffset = 0,
		leftFunc = 0, rightFunc = 0,
        topBlendLeftStart = 0, topBlendRightStart = 0,
        topBlendPeakOffset = 0, topBlendPeak = 0,
        topBlendOffset = 0, topBlendLeftOffset = 0, topBlendRightOffset = 0,
        xIt, yIt;
	float blendLeftA = 0, blendRightA = 0, curEU = 0, topBlendA = 0,
        leftA, rightA;
	coord startFill = coord(), endFill = coord(),
		a = coord(), b = coord(), c = coord(), d = coord();

	// 0 = exponent decay, 1 = linear + blend, 2 = parabola + blend
    /*parabola n exp based on skew*/
    /*if (skew == 0.0) leftFunc = 1, rightFunc = 1,
        leftA = LinearM(curCoord.height, radius * 0.75),
        rightA = leftA,
        blendLeftA = ExponentDecayA(curCoord.height, radius * 0.75),
        blendRightA = blendLeftA,
        blendLeftOffset = LinearX(leftA, curCoord.height, curCoord.height / 4),
        blendRightOffset = blendLeftOffset,
        blendLeftHeight = (curCoord.height / 4) / blendLeftA,
        blendRightHeight = blendLeftHeight;
	else if (skew > 0.0) leftFunc = 0, rightFunc = 2,
		leftA = ExponentDecayA(curCoord.height, leftRadius),
		rightA = ParabolaA(curCoord.height, rightRadius * 0.75),
		blendRightA = ExponentDecayA(curCoord.height, rightRadius * 0.75),
		blendRightOffset = ParabolaX(rightA, curCoord.height, curCoord.height / 4),
        blendRightHeight = (curCoord.height / 4) / blendRightA;
	else leftFunc = 2, rightFunc = 0,
		leftA = ParabolaA(curCoord.height, leftRadius * 0.75),
		blendLeftA = ExponentDecayA(curCoord.height, leftRadius * 0.75),
		blendLeftOffset = ParabolaX(leftA, curCoord.height, curCoord.height / 4),
        blendLeftHeight = (curCoord.height / 4) / blendLeftA,
		rightA = ExponentDecayA(curCoord.height, rightRadius);*/
    /*always linear (+ exponent blend)*/
    leftA = LinearM(curCoord.height, leftRadius * 0.75),
        blendLeftA = ExponentDecayA(curCoord.height, leftRadius * 0.75),
        blendLeftOffset = LinearX(leftA, curCoord.height, curCoord.height / 4),
        blendLeftHeight = (curCoord.height / 4) / blendLeftA;
    rightA = LinearM(curCoord.height, rightRadius * 0.75),
        blendRightA = ExponentDecayA(curCoord.height, rightRadius * 0.75),
        blendRightOffset = LinearX(rightA, curCoord.height, curCoord.height / 4),
        blendRightHeight = (curCoord.height / 4) / blendRightA;
    topBlendLeftStart = topBlend * leftRadius,
        topBlendRightStart = topBlend * rightRadius,
        topBlendPeakOffset = Linear(leftA, topBlendLeftStart, curCoord.height),
        topBlendPeak = (curCoord.height - topBlendPeakOffset) / 2,
        topBlendOffset = (topBlendLeftStart + topBlendRightStart) / 2,
        topBlendA = ParabolaA(topBlendPeak, topBlendOffset),
        topBlendLeftOffset = topBlendLeftStart - topBlendOffset,
        topBlendRightOffset = topBlendRightStart - topBlendOffset;

	TArray<coord> grad;

	// figuring out where to start loop;
	a.x = 
        FMath::RoundToInt(FMath::Sin(DegreeToRad(curCoord.AddThetaTemp(180))) * radius) +
        FMath::RoundToInt(FMath::Sin(DegreeToRad(curCoord.AddThetaTemp(-90))) * leftRadius);
	a.y = 
        FMath::RoundToInt(FMath::Cos(DegreeToRad(curCoord.AddThetaTemp(180))) * radius) +
        FMath::RoundToInt((float)FMath::Cos(DegreeToRad(curCoord.AddThetaTemp(-90))) * leftRadius);
	b.x = 
        FMath::RoundToInt(FMath::Sin(DegreeToRad(curCoord.AddThetaTemp(180))) * radius) +
        FMath::RoundToInt(FMath::Sin(DegreeToRad(curCoord.AddThetaTemp(90))) * rightRadius);
	b.y = 
        FMath::RoundToInt(FMath::Cos(DegreeToRad(curCoord.AddThetaTemp(180))) * radius) +
        FMath::RoundToInt(FMath::Cos(DegreeToRad(curCoord.AddThetaTemp(90))) * rightRadius);
	c.x = 
        FMath::RoundToInt(FMath::Sin(DegreeToRad(curCoord.theta)) * radius) +
        FMath::RoundToInt(FMath::Sin(DegreeToRad(curCoord.AddThetaTemp(90))) * rightRadius);
	c.y = 
        FMath::RoundToInt(FMath::Cos(DegreeToRad(curCoord.theta)) * radius) +
        FMath::RoundToInt(FMath::Cos(DegreeToRad(curCoord.AddThetaTemp(90))) * rightRadius);
	d.x = 
        FMath::RoundToInt(FMath::Sin(DegreeToRad(curCoord.theta)) * radius) +
        FMath::RoundToInt(FMath::Sin(DegreeToRad(curCoord.AddThetaTemp(-90))) * leftRadius);
	d.y = 
        FMath::RoundToInt(FMath::Cos(DegreeToRad(curCoord.theta)) * radius) +
        FMath::RoundToInt(FMath::Cos(DegreeToRad(curCoord.AddThetaTemp(-90))) * leftRadius);

	// compare Abs value; multply by left n right radius in order for skew to work properly
	startFill.x = (Abs(a.x) > Abs(d.x)) ? a.x : d.x;
	startFill.y = (Abs(a.y * rightRadius) > Abs(b.y * leftRadius)) ? a.y : b.y;
	endFill.x = (Abs(c.x) > Abs(b.x)) ? c.x : b.x;
	endFill.y = (Abs(c.y * leftRadius) > Abs(d.y * rightRadius)) ? c.y : d.y;

    xIt = (endFill.x > startFill.x) ? 1 : -1;
	yIt = (endFill.y > startFill.y) ? 1 : -1;

    grad.Init(coord(), Abs((endFill.x - startFill.x + xIt) * (endFill.y - startFill.y + yIt)));

    for (int x = startFill.x; (xIt == 1) ? x <= endFill.x : x >= endFill.x; x += xIt) {
        for (int y = startFill.y; (yIt == 1) ? y <= endFill.y : y >= endFill.y; y += yIt) {
            curIndex = Abs((x - startFill.x) * (endFill.y - startFill.y)) + Abs(y - startFill.y);
            // fill in x degree left and right only
            tempTheta = RadToDegree(FMath::Atan2(y, x));
            tempTheta = -(tempTheta - 90);
            grad[curIndex].SetTheta(tempTheta); // 0 @ north
            grad[curIndex].AddTheta(-curCoord.theta); // 0 @ center heading
            if ((grad[curIndex].theta > 270 - fillDegree && grad[curIndex].theta < 270 + fillDegree) || (x == 0 && y == 0)) {
                curEU = EuclideanDistance(coord(x, y, 0)); // calculate if condition match only; optimization
                // left
                grad[curIndex].x = x + curCoord.x;
                grad[curIndex].y = y + curCoord.y;

                //grad[curIndex].height = 20;
                /*switch (leftFunc) {
                case 0:
                    grad[curIndex].height = ExponentDecay(leftA, curCoord.height, curEU);
                    break;
                case 1:
                    if (curEU > blendLeftOffset) grad[curIndex].height = ExponentDecay(blendLeftA, blendLeftHeight, curEU, blendLeftOffset);
                    else grad[curIndex].height = Linear(leftA, curEU, curCoord.height);
                    break;
                case 2:
                    if (curEU > blendLeftOffset) grad[curIndex].height = ExponentDecay(blendLeftA, blendLeftHeight, curEU, blendLeftOffset);
                    else grad[curIndex].height = Parabola(leftA, curCoord.height, curEU);
                    break;
                }*/
                /*linear only*/
                if (curEU > blendLeftOffset) grad[curIndex].height = ExponentDecay(blendLeftA, blendLeftHeight, curEU, blendLeftOffset);
                else if (curEU < topBlendLeftStart) grad[curIndex].height = Parabola(topBlendA, topBlendPeak + topBlendPeakOffset, Abs(curEU + topBlendLeftOffset));
                else grad[curIndex].height = Linear(leftA, curEU, curCoord.height);
            }
            else if (grad[curIndex].theta > 90 - fillDegree && grad[curIndex].theta < 90 + fillDegree) {
                curEU = EuclideanDistance(coord(x, y, 0)); // calculate if condition match only; optimization
                // right
                grad[curIndex].x = x + curCoord.x;
                grad[curIndex].y = y + curCoord.y;
                //grad[curIndex].height = 20;
                /*switch (rightFunc) {
                case 0:
                    grad[curIndex].height = ExponentDecay(rightA, curCoord.height, curEU);
                    break;
                case 1:
                    if (curEU > blendRightOffset) grad[curIndex].height = ExponentDecay(blendRightA, blendRightHeight, curEU, blendRightOffset);
                    else grad[curIndex].height = Linear(rightA, curEU, curCoord.height);
                    break;
                case 2:
                    if (curEU > blendRightOffset) grad[curIndex].height = ExponentDecay(blendRightA, blendRightHeight, curEU, blendRightOffset);
                    else grad[curIndex].height = Parabola(rightA, curCoord.height, curEU);
                    break;
                }*/
                /*linear only*/
                if (curEU > blendRightOffset) grad[curIndex].height += ExponentDecay(blendRightA, blendRightHeight, curEU, blendRightOffset);
                else if (curEU < topBlendRightStart) grad[curIndex].height += Parabola(topBlendA, topBlendPeak + topBlendPeakOffset, Abs(curEU + topBlendRightOffset));
                else grad[curIndex].height += Linear(rightA, curEU, curCoord.height);
            }
        }
    }
    isAdd ? DrawDetail(grad) : Draw(grad);
}

float ULanGenElevation::EuclideanDistance(coord pointCoord, coord centerCoord)
{
    return FMath::Sqrt(
        FMath::Pow(centerCoord.x - pointCoord.x, 2) + FMath::Pow(centerCoord.y - pointCoord.y, 2)
    );
}

int ULanGenElevation::Parabola(float a, float c, float x, float xOffset)
{
    x -= xOffset;
    return (a * x * x) + c;
}

float ULanGenElevation::ParabolaA(float c, float xMax, float xOffset)
{
    xMax -= xOffset;
    return (float)-c / (xMax * xMax);
}

float ULanGenElevation::ParabolaX(float a, float c, float y)
{
    return FMath::Pow((y - c)/a , 0.5);
}

int ULanGenElevation::ExponentDecay(float a, float b, float x, float xOffset, int modifier)
{
    x -= xOffset;
    return b * FMath::Pow(a, x * modifier);
}

float ULanGenElevation::ExponentDecayA(float b, float xMax, float y, float xOffset)
{
    xMax -= xOffset;
    return FMath::Pow(y / b, 1 / xMax);
}

int ULanGenElevation::Linear(float m, float x, float c) { return m * x + c; }

float ULanGenElevation::LinearM(float c, float xMax, float modifier) { return modifier * (c / xMax); }

float ULanGenElevation::LinearX(float m, float c, float y) { return (y - c) / m; }

TArray<coord> ULanGenElevation::BezierCurve3(coord p1, coord p2, coord p3)
{
    TArray<coord> res;
    float step = 0.1, x, y;
    //P = (1−t)2P1 + 2(1−t)tP2 + t2P3
    for (float t = 0; t <= 1.0; t += step) {
        x = ((1 - t) * (1 - t)) * p1.x +
            2 * (1 - t) * t * p2.x +
            (t * t) * p3.x;
        y = ((1 - t) * (1 - t)) * p1.y +
            2 * (1 - t) * t * p2.y +
            (t * t) * p3.y;
        res.Add(coord(x, y, 0));
    }
    return res;
}

TArray<coord> ULanGenElevation::BezierCurve4(coord p1, coord p2, coord p3, coord p4)
{
    TArray<coord> res;
    float step = 0.1, x, y;
    //P = (1−t)3P1 + 3(1−t)2tP2 + 3(1−t)t2P3 + t3P4
    for (float t = 0; t <= 1.0; t += step) {
        x = ((1 - t) * (1 - t) * (1 - t)) * p1.x +
            3 * ((1 - t) * (1 - t)) * t * p2.x +
            3 * (1 - t) * (t * t) * p3.x +
            (t * t * t) * p4.x;
        y = ((1 - t) * (1 - t) * (1 - t)) * p1.y +
            3 * ((1 - t) * (1 - t)) * t * p2.y +
            3 * (1 - t) * (t * t) * p3.y +
            (t * t * t) * p4.y;
        res.Add(coord(x, y, 0));
    }
    return res;
}

void ULanGenElevation::Draw(TArray<coord> currentLine, bool overwrite)
{
    /*0 = do not overwerite; 1 = overwrite; 2 = add value*/
    for (coord i : currentLine) {
        if (i.isInRange(lanX, lanY)) {
            /*only draw if current height higher*/
            if (i.height > texture[i.index(lanY)].R - init.R || overwrite) texture[i.index(lanY)].R = i.height + init.R;
        }
    }
}

void ULanGenElevation::DrawDetail(TArray<coord> currentLine, bool overwrite)
{
    /*0 = do not overwerite; 1 = overwrite; 2 = add value*/
    for (coord i : currentLine) {
        if (i.isInRange(lanX, lanY)) {
            /*only draw if current height higher*/
            if (i.height > detailTexture[i.index(lanY)].R || overwrite) detailTexture[i.index(lanY)].R = i.height;
        }
    }
}

float ULanGenElevation::DegreeToRad(int degree) { return degree * 3.14159265 / 180; }

float ULanGenElevation::RadToDegree(float rad) { return rad * 180 / 3.14159265; }

int ULanGenElevation::Abs(int in) { return (in < 0) ? in * -1 : in; }