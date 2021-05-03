/*
* Bresenham converted from https://gist.github.com/polyamide/424a9fc8507436883a911e6f1d79d9d3#file-bresenham-js-L1
*
* Edited and converted to Unreal compatible c++ by Fachrurrozy Muhammad
*/

#include "LanGenElevation.h"
#include "Math/UnrealMathUtility.h"
#include "Misc/Char.h"
#include "Misc/DefaultValueHelper.h"

#define PI 3.14159265
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

TArray<FColor> ULanGenElevation::GenerateGraph(
    FString rule, FString axiom, int ruleLoop, 
    int midPointLoop, int midPointSmooth, int x,
    int y, int lineLength, int minAngle,
    int maxAngle, bool generateElevation, int radius,
    int peak, int gapAccuracy)
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
        case 'F' || 'G' || 'H': Bresenham(currentLine, lineLength); break;
        case 'P': /*peak point*/ break;
        case 'L': /*lowest point*/ break;
        case '+': currentCoord->AddTheta(isRandomAngle ? randomEngine.RandRange(minAngle, maxAngle) : minAngle); break;
        case '-': currentCoord->AddTheta((isRandomAngle ? randomEngine.RandRange(minAngle, maxAngle) : minAngle) * -1); break;
        case '[': branchRootStack.Add(*currentCoord); break;
        case ']': // run on line ends
            MidpointDisplacement(currentLine, peak, midPointSmooth, midPointLoop);
            //Draw(currentLine);
            //GradientFill(currentLine, radius, peak, gapAccuracy);
            for (coord j : currentLine) CiircularGradientFill(j, radius);
            currentLine.Empty();
            currentLine.Add(
                branchRootStack[branchRootStack.Num() - 1]
            );
            branchRootStack.RemoveAt(branchRootStack.Num() - 1);
            break;
        }
    }
    MidpointDisplacement(currentLine, peak, midPointSmooth, midPointLoop);
    //Draw(currentLine);
    //GradientFill(currentLine, radius, peak, gapAccuracy);
    for (coord i : currentLine) CiircularGradientFill(i, radius);

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
    int x = x0 + FMath::Sin(currentCoord.theta * PI / 180) * lineLength;
    int y = y0 + FMath::Cos(currentCoord.theta * PI / 180) * lineLength;

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

void ULanGenElevation::GradientFill(TArray<coord>& currentLine, int radius, int peak, int gapAccuracy)
{
    // radius @ peak
    //CON_LOG("GradientFill ver 42");
    float currentRadius = radius;
    float a = FMath::Pow((float)radius, (1 / (float)peak)); // exponent ridge height
    float b = currentRadius / peak; // linear ridge height
    //CON_LOG("currLine size, a, peak, radius: %d, %f, %d, %d", currentLine.Num() , a, peak, radius);
    TArray<coord> left, right;
    coord currentCoord;
    for (int i = 0; i < currentLine.Num(); ++i) {
        // exponent
        //currentRadius = FMath::Pow(a, currentLine[i].height);
        // linear
        currentRadius = b * currentLine[i].height;

        // check for angle change
        if (i > 0) {
            if (currentLine[i].theta != currentLine[i - 1].theta) {
                CON_LOG("different angle");
                // do angle gap fill
                if (currentLine[i].theta > currentLine[i - 1].AddThetaTemp(180) || currentLine[i].theta < currentLine[i - 1].theta) {
                    // go left, fill right
                    CON_LOG("go left | current: %d, before: %d version 2", currentLine[i].theta, currentLine[i-1].theta);
                    int j = currentLine[i].theta;
                    while (j < currentLine[i - 1].theta) {
                        currentCoord = currentLine[i];
                        currentCoord.theta = j;
                        right.Empty();
                        currentCoord.AddTheta(90);
                        right.Add(currentCoord);
                        Bresenham(right, currentRadius);
                        // gradient (linear)
                        right[right.Num() - 1].height = 0;
                        for (int k = 0; k < right.Num() - 1; ++k) {
                            right[k].height = ((float)k / (right.Num() - 1)) * ((float)right[right.Num() - 1].height - right[0].height) + right[0].height;
                            //if (i == 50) CON_LOG("batch %d height: %d rad: %f size: %d", i, right[k].height, currentRadius, right.Num());
                        }
                        Draw(right);
                        j += gapAccuracy;
                    }
                }
                else {
                    // go right, fill left
                    CON_LOG("go right | current: %d, before: %d", currentLine[i].theta, currentLine[i - 1].theta);
                    int j = currentLine[i].theta;
                    while (j > currentLine[i - 1].theta) {
                        currentCoord = currentLine[i];
                        currentCoord.theta = j;
                        left.Empty();
                        currentCoord.AddTheta(-90);
                        left.Add(currentCoord);
                        Bresenham(left, currentRadius);
                        // gradient (linear)
                        left[left.Num() - 1].height = 0;
                        for (int k = 0; k < left.Num() - 1; ++k) {
                            left[k].height = ((float)k / (left.Num() - 1)) * ((float)left[left.Num() - 1].height - left[0].height) + left[0].height;
                            //if (i == 50) CON_LOG("batch %d height: %d rad: %f size: %d", i, left[k].height, currentRadius, left.Num());
                        }
                        Draw(left);
                        j -= gapAccuracy;
                    }
                }
            }
        }

        // do radius fill | L R balanced radius feature ?
        currentCoord = currentLine[i];

        // left side
        left.Empty();
        currentCoord.AddTheta(-90);
        left.Add(currentCoord);
        //CON_LOG("pre-Bresenham, left size: %d, %d", currentCoord.height, left.Num());
        Bresenham(left, currentRadius);
        // gradient (linear)
        left[left.Num() - 1].height = 0;
        for (int k = 0; k < left.Num() - 1; ++k) {
            left[k].height = ((float)k / (left.Num() - 1)) * ((float)left[left.Num() - 1].height - left[0].height) + left[0].height;
            //if (i == 50) CON_LOG("batch %d height: %d rad: %f size: %d", i, left[k].height, currentRadius, left.Num());
        }
        Draw(left);

        // right side
        right.Empty();
        currentCoord.AddTheta(180);
        right.Add(currentCoord);
        //CON_LOG("pre-Bresenham, right size: %d, %d", currentCoord.height, right.Num());
        Bresenham(right, currentRadius);
        //// gradient (linear)
        right[right.Num() - 1].height = 0;
        for (int k = 0; k < right.Num() - 1; ++k) {
            right[k].height = ((float)k / (right.Num() - 1)) * ((float)right[right.Num() - 1].height - right[0].height) + right[0].height;
        }
        Draw(right);
    }
}

void ULanGenElevation::CiircularGradientFill(coord centerCoord, int radius)
{
    TArray<coord> gradient;
    TArray<float> unnormalized;
    int size = radius * 2 + 1,
        counter = 0,
        peak = centerCoord.height;
    float biggest = 0;
    coord startCoord = coord(centerCoord.x - radius, centerCoord.y - radius, 0);
    //CON_LOG("center: %d, %d; start: %d, %d", centerCoord.x, centerCoord.y, startCoord.x, startCoord.y);

    gradient.Init(coord(), size * size);
    unnormalized.Init(0, size * size);
    // coord sync
    for (coord& i : gradient) {
        i.x = (int)counter / size + startCoord.x, i.y = (int)counter % size + startCoord.y, ++counter;
        //CON_LOG("gradient: %d, %d", i.x, i.y);
    }
    // gradient value
    for (int row = 0; row < size; ++row) {
        for (int col = 0; col < size; ++col) {
            unnormalized[row * size + col] = EuclideanDistance(centerCoord, gradient[row * size + col], radius);
            if (unnormalized[row * size + col] > biggest) biggest = unnormalized[row * size + col];
        }
    }
    //CON_LOG("biggest: %f", biggest);
    // normalized
    for (int i = 0; i < gradient.Num(); ++i) {
        // invert then normalize
        if (unnormalized[i] > 0) gradient[i].height = ((float)((biggest - unnormalized[i]) / biggest) * peak);
    }

    Draw(gradient);
}

void ULanGenElevation::CiircularLineGradientFill(TArray<coord> currentLine, int radius)
{
    
}

float ULanGenElevation::EuclideanDistance(coord centerCoord, coord pointCoord, int radius)
{
    float dist = FMath::Sqrt(
        FMath::Pow(centerCoord.x - pointCoord.x, 2) + FMath::Pow(centerCoord.y - pointCoord.y, 2)
    );
    return (dist > radius ? 0 : dist);
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
