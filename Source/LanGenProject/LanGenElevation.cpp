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
    int x, int y, int lineLength, int minAngle, 
    int maxAngle, bool generateElevation, int radius,
    int peak)
{
    lanX = x, lanY = y;
    TArray<FColor> texture;
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
        case '-': currentCoord->AddTheta((isRandomAngle ? randomEngine.RandRange(minAngle, maxAngle) : minAngle) * -1); break;
        case '[': branchRootStack.Add(*currentCoord); break;
        case ']': // run on line ends
            MidpointDisplacement(currentLine, peak);
            Draw(texture, currentLine);
            GradientFill(texture, currentLine, radius, peak);
            currentLine.Empty();
            currentLine.Add(
                branchRootStack[branchRootStack.Num() - 1]
            );
            /*CON_LOG("branch start: %d %d %d", 
                branchRootStack[branchRootStack.Num() - 1].x,
                branchRootStack[branchRootStack.Num() - 1].y,
                branchRootStack[branchRootStack.Num() - 1].theta);*/
            branchRootStack.RemoveAt(branchRootStack.Num() - 1);
            break;
        }
    }
    MidpointDisplacement(currentLine);
    Draw(texture, currentLine);
    GradientFill(texture, currentLine, radius, peak);

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
    int x = x0 + sin(currentCoord.theta * PI / 180) * lineLength;
    int y = y0 + cos(currentCoord.theta * PI / 180) * lineLength;

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
void ULanGenElevation::MidpointDisplacement(TArray<coord>& currentLineCoord, int loop, int displacement)
{
    TArray<midPoint> oldIndexes, newIndexes;
    midPoint mid;
    int currentDisplacement = displacement,
        randomModifier;
    float smooth = 0.45;

    if (FMath::Pow(2, loop) >= currentLineCoord.Num()) loop = FMath::Log2(currentLineCoord.Num());
    oldIndexes.Add(midPoint(0, currentLineCoord[0].height));
    oldIndexes.Add(midPoint(currentLineCoord.Num() - 1, currentLineCoord[currentLineCoord.Num() - 1].height));

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
                currentDisplacement *= smooth;
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

void ULanGenElevation::GradientFill(TArray<FColor>& texture, TArray<coord> currentLine, int radius, int peak)
{
    // radius @ peak
    CON_LOG("Gradient fill called");
    float currentRadius = radius,
        a = FMath::Pow((float)radius, (1 / (float)peak));
    TArray<coord> left, right;
    coord currentCoord;
    for (int i = 0; i < currentLine.Num(); ++i) {
        currentRadius = FMath::Pow(a, currentLine[i].height);
        CON_LOG("radius pow: %f", currentRadius);
        if (i != 0) {
            // check for angle change
            if (currentLine[i].theta != currentLine[i - 1].theta) {
                CON_LOG("different angle");
                // do angle gap fill
                if (currentLine[i].theta > currentLine[i - 1].theta + 180 || currentLine[i].theta < currentLine[i - 1].theta) {
                    // go left, fill right
                    for (int j = currentLine[i].theta;
                        j > currentLine[i - 1].theta;
                        j += (currentLine[i].theta > currentLine[i - 1].theta ? 10 : -10)) {
                        currentCoord = currentLine[j];
                        right.Empty();
                        currentCoord.AddTheta(90);
                        right.Add(currentCoord);
                        Bresenham(right, currentRadius);

                        // gradient (linear)
                        right[right.Num() - 1].height = 0;
                        for (int k = 0; k < right.Num() - 1; ++k) {
                            right[k].height = ((float) k / (right.Num() - 1)) * (right[right.Num() - 1].height - right[0].height);
                        }

                        Draw(texture, right);
                    }
                }
                else {
                    // go right, fill left
                    for (int j = currentLine[i].theta;
                        j > currentLine[i - 1].theta;
                        i += (currentLine[i].theta > currentLine[i - 1].theta ? 10 : -10)) {
                        currentCoord = currentLine[j];
                        left.Empty();
                        currentCoord.AddTheta(-90);
                        left.Add(currentCoord);
                        Bresenham(left, currentRadius);

                        // gradient (linear)
                        left[left.Num() - 1].height = 0;
                        for (int k = 0; k < left.Num() - 1; ++k) {
                            left[k].height = ((float)k / (left.Num() - 1)) * (left[left.Num() - 1].height - left[0].height);
                        }

                        Draw(texture, left);
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
        Bresenham(left, currentRadius);

        // gradient (linear)
        left[left.Num() - 1].height = 0;
        /*CON_LOG("left size: %d", left.Num());*/
        for (int k = 0; k < left.Num() - 1; ++k) {
            left[k].height = ((float)k / (left.Num() - 1)) * (left[left.Num() - 1].height - left[0].height);
        }

        Draw(texture, left);

        // right side
        right.Empty();
        currentCoord.AddTheta(180);
        right.Add(currentCoord);
        Bresenham(right, currentRadius);

        // gradient (linear)
        right[right.Num() - 1].height = 0;
        for (int k = 0; k < right.Num() - 1; ++k) {
            right[k].height = ((float)k / (right.Num() - 1)) * (right[right.Num() - 1].height - right[0].height);
        }

        Draw(texture, right);
    }
}

void ULanGenElevation::Draw(TArray<FColor>& texture, TArray<coord> currentLine)
{
    //CON_LOG("height: %d", currentLine[0].height);
    for (coord i : currentLine) {
        if (i.x >= 0 && i.y >= 0 && i.x < lanX && i.y < lanY) 
            /*only draw if current height higher*/
            if (i.height > texture[i.x * lanY + i.y].R - init.R) texture[i.x * lanY + i.y].R = i.height + init.R;
    }
}
