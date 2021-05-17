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

void ULanGenElevation::InitSeed(int32 in)
{
	seed = in;
	p.Empty();
	p = P_BASE;
	randomEngine = FRandomStream(seed);
	TArray<int> tempArray = P_BASE;
	Shuffle(tempArray);
	p = tempArray;
	p.Append(tempArray);
}

TArray<FColor> ULanGenElevation::GenerateGraphDebug(
    FString rule, FString axiom, int ruleLoop,
    int x, int y, int lineLength,
    int minAngle, int maxAngle, int radius,
    int peak, float skew, int fillDegree
)
{
    lanX = x, lanY = y;
    TArray<coord> branchRootStack;
    TArray<coord> currentLine;
    coord* currentCoord;
    FString grammar;
    bool isRandomAngle = minAngle != maxAngle;

    init = FColor(100, 0, 0);
    texture.Init(init, lanX * lanY);

    // L-System
    RuleSetup(rule);
    grammar = RuleApply(axiom, ruleLoop);
    //CON_LOG("%s", *grammar);

    // set starting position
    currentLine.Add(coord(
        lanX / 2,
        lanY / 2,
        0
    ));

    // create array of target coord
    for (TCHAR i : grammar) {
        currentCoord = &currentLine[currentLine.Num() - 1];
        switch (i) {
        case 'F': Bresenham(currentLine, lineLength); break;
        case 'P': /*peak point*/ break;
        case 'L': /*lowest point*/ break;
        case '+': currentCoord->AddTheta(isRandomAngle ? randomEngine.RandRange(minAngle, maxAngle) : minAngle); break;
        case '-': currentCoord->AddTheta(-(isRandomAngle ? randomEngine.RandRange(minAngle, maxAngle) : minAngle)); break;
        case '[': branchRootStack.Add(*currentCoord); break;
        case ']': // run on line ends
            GradientSingleMain(currentLine, peak, radius, skew, fillDegree);
            currentLine.Empty();
            currentLine.Add(
                branchRootStack[branchRootStack.Num() - 1]
            );
            branchRootStack.RemoveAt(branchRootStack.Num() - 1);
            break;
        }
    }
    GradientSingleMain(currentLine, peak, radius, skew, fillDegree);

    return texture;

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
    FString currentLstring;
    for (int i = 0; i < loop; ++i) {
        currentLstring = "";
        //check for each char in axiom
        for (TCHAR j : axiom) {
            // rule check
            for (int k = 0; k < rules.Num(); ++k) {
                if (j == rules[k].from[0]) {
                    currentLstring.Append(RandomizeRule(k));
                    break;
                }
            }
            currentLstring.AppendChar(j);
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
FString ULanGenElevation::RandomizeRule(int ruleIndex)
{
    int randomNumber = randomEngine.RandRange(1, 100);
    int currentRequirement = 0;
    for (int i = 0; i < rules[ruleIndex].to.Num(); ++i) {
        currentRequirement += rules[ruleIndex].prob[i];
        if (randomNumber <= currentRequirement) {
            //CON_LOG("%s", *rules[ruleIndex].to[i]);
            return rules[ruleIndex].to[i];
        }
    }
    return rules[ruleIndex].from;
}

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
void ULanGenElevation::MidpointDisplacement(TArray<coord>& currentLineCoord, int displacement, int smooth, int loop)
{
    CON_LOG("midpoint");
    TArray<midPoint> oldIndexes, newIndexes;
    midPoint mid;
    int currentDisplacement = displacement,
        randomModifier;
    float modifier = FMath::Pow(2, -smooth);

    if (FMath::Pow(2, loop) >= currentLineCoord.Num()) {
        loop = FMath::Log2(currentLineCoord.Num());
        //CON_LOG("loop chopped");
    }
    //CON_LOG("current loop: %d loop result: %d limit: %d", loop, (int)FMath::Pow(2, loop), currentLineCoord.Num());
    // if within landscape
    if (currentLineCoord[0].isInRange(lanX, lanY)) {
        oldIndexes.Add(midPoint(0,
            texture[currentLineCoord[0].index(lanY)].R - init.R));
    }
    else oldIndexes.Add(midPoint(0, 0));

    if (currentLineCoord[currentLineCoord.Num() - 1].isInRange(lanX, lanY)) {
        oldIndexes.Add(midPoint(currentLineCoord.Num() - 1,
            texture[currentLineCoord[currentLineCoord.Num() - 1].index(lanY)].R - init.R));
    }
    else oldIndexes.Add(midPoint(currentLineCoord.Num() - 1, 0));

    // if branch
    if (oldIndexes[0].height > 0) currentDisplacement *= modifier;

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

    for (int i = 0; i < loop; ++i) {
        // mid point
        for (int j = 0; j < oldIndexes.Num() - 1; ++j) {
            mid.index = (int)((oldIndexes[j].index + oldIndexes[j + 1].index) / 2);
            newIndexes.Add(oldIndexes[j]);
            if (mid.index != oldIndexes[j].index) {
                if (currentLineCoord[mid.index].height < currentDisplacement) randomModifier = 1;
                else randomModifier = (randomEngine.RandRange(0, 1) == 0 ? 1 : -1);
                mid.height = currentLineCoord[mid.index].height + currentDisplacement * randomModifier;
                currentDisplacement *= modifier;
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

    /*CON_LOG("Line Start =================");
    for (coord& i : currentLineCoord) {
        CON_LOG("height: %d", i.height);
    }*/
}

void ULanGenElevation::GradientSingleMain(TArray<coord>& curLine, int peak, int radius, float skew, int fillDegree)
{
    //CON_LOG("grad single main called v8");

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
    // gradient; peak -> x; radius -> y; y = mx; m = y / x
    int eachRadius;
    float m = (float) radius / peak;

    Draw(curLine);

    for (coord curCoord : curLine) {
        eachRadius = m * curCoord.height;
        //CON_LOG("m calc: %f \teach rad: %d, \theight: %d", (float) radius / peak, eachRadius, curCoord.height);
        GradientSingleMainHelper(curCoord, eachRadius, skew, fillDegree);
    }
}

void ULanGenElevation::GradientSingleMainHelper(coord curCoord, int radius, float skew, int fillDegree)
{
	int leftRadius = radius * ((-skew) + 1),
		rightRadius = radius * (skew + 1),
        curIndex = 0, tempTheta = 0,
		leftFunc, rightFunc, xIt, yIt, blendOffset, blendHeight;
	float leftA, rightA, optiA, curEU = 0;
	coord startFill = coord(), endFill = coord(),
		a = coord(), b = coord(), c = coord(), d = coord();

	// 0 = exponent decay, 1 = linear + blend, 2 = parabola + blend
	if (skew == 0.0) leftFunc = 1, rightFunc = 1,
		leftA = LinearM(curCoord.height, radius * 0.75),
		rightA = leftA,
		optiA = ExponentDecayA(curCoord.height, radius * 0.75),
		blendOffset = LinearX(leftA, curCoord.height, curCoord.height / 4);
	else if (skew > 0.0) leftFunc = 0, rightFunc = 2,
		leftA = ExponentDecayA(curCoord.height, leftRadius),
		rightA = ParabolaA(curCoord.height, rightRadius * 0.75),
		optiA = ExponentDecayA(curCoord.height, rightRadius * 0.75),
		blendOffset = ParabolaX(rightA, curCoord.height, curCoord.height / 4);
	else leftFunc = 2, rightFunc = 0,
		leftA = ParabolaA(curCoord.height, leftRadius * 0.75),
		optiA = ExponentDecayA(curCoord.height, leftRadius * 0.75),
		blendOffset = ParabolaX(leftA, curCoord.height, curCoord.height / 4),
		rightA = ExponentDecayA(curCoord.height, rightRadius);

	blendHeight = (curCoord.height / 4) / optiA;

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

    /*CON_LOG("OwO start (%d, %d) \tend (%d, %d) \tsize %d", startFill.x, startFill.y, endFill.x, endFill.y,
        Abs((endFill.x - startFill.x + xIt) * (endFill.y - startFill.y + yIt))
    );*/
    /*CON_LOG("a (%d, %d); b (%d, %d); c (%d, %d); d (%d, %d); itt (%d, %d); rad (%d, %d, %d)", a.x, a.y, b.x, b.y, c.x, c.y, d.x, d.y, xIt, yIt,
        leftRadius, radius, rightRadius
    );*/

	for (int x = startFill.x; (xIt == 1) ? x <= endFill.x : x >= endFill.x; x += xIt) {
		for (int y = startFill.y; (yIt == 1) ? y <= endFill.y : y >= endFill.y; y += yIt) {
			curIndex = Abs((x - startFill.x) * (endFill.y - startFill.y)) + Abs(y - startFill.y);
			// fill in x degree left and right only
            tempTheta = RadToDegree(FMath::Atan2(y, x));
            tempTheta = -(tempTheta - 90);
            grad[curIndex].SetTheta(tempTheta); // 0 @ north
            grad[curIndex].AddTheta(-curCoord.theta); // 0 @ center heading
			
			if (grad[curIndex].theta > 270 - fillDegree && grad[curIndex].theta < 270 + fillDegree) {
                curEU = EuclideanDistance(coord(x, y, 0)); // calculate if condition match only; optimization
				// left
				grad[curIndex].x = x + curCoord.x;
				grad[curIndex].y = y + curCoord.y;
                //grad[curIndex].height = 20;
				switch (leftFunc) {
				case 0:
					grad[curIndex].height = ExponentDecay(leftA, curCoord.height, curEU);
					break;
				case 1:
					if (curEU > blendOffset) grad[curIndex].height = ExponentDecay(optiA, blendHeight, curEU, blendOffset);
					else grad[curIndex].height = Linear(leftA, curEU, curCoord.height);
					break;
				case 2:
					if (curEU > blendOffset) grad[curIndex].height = ExponentDecay(optiA, blendHeight, curEU, blendOffset);
					else grad[curIndex].height = Parabola(leftA, curCoord.height, curEU);
					break;
				}
			}
			else if (grad[curIndex].theta > 90 - fillDegree && grad[curIndex].theta < 90 + fillDegree) {
                curEU = EuclideanDistance(coord(x, y, 0));
				// right
				grad[curIndex].x = x + curCoord.x;
				grad[curIndex].y = y + curCoord.y;
                //grad[curIndex].height = 20;
				switch (rightFunc) {
				case 0:
					grad[curIndex].height = ExponentDecay(rightA, curCoord.height, curEU);
					break;
				case 1:
					if (curEU > blendOffset) grad[curIndex].height = ExponentDecay(optiA, blendHeight, curEU, blendOffset);
					else grad[curIndex].height = Linear(rightA, curEU, curCoord.height);
					break;
				case 2:
					if (curEU > blendOffset) grad[curIndex].height = ExponentDecay(optiA, blendHeight, curEU, blendOffset);
					else grad[curIndex].height = Parabola(rightA, curCoord.height, curEU);
					break;
				}
			}
		}
	}
	Draw(grad);
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

int ULanGenElevation::Linear(float m, float x, float c)
{
    return m * x + c;
}

float ULanGenElevation::LinearM(float c, float xMax)
{
    return -(FMath::Pow(2, (c / xMax) - 1));
}

float ULanGenElevation::LinearX(float m, float c, float y)
{
    return (y - c) / m;
}

void ULanGenElevation::Draw(TArray<coord> currentLine)
{
    //CON_LOG("height: %d", currentLine[0].height);
    for (coord i : currentLine) {
        if (i.isInRange(lanX, lanY)) 
            /*only draw if current height higher*/
            if (i.height > texture[i.index(lanY)].R - init.R) texture[i.index(lanY)].R = i.height + init.R;
    }
}

float ULanGenElevation::DegreeToRad(int degree)
{
    return degree * 3.14159265 / 180;
}

float ULanGenElevation::RadToDegree(float rad)
{
    return rad * 180 / 3.14159265;
}

int ULanGenElevation::Abs(int in)
{
    return (in < 0) ? in * -1 : in;
}