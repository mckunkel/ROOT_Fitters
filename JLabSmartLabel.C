void JLabSmartLabel(Double_t xpos=0.9, Double_t ypos=0.9, Double_t scale=1.0, TString str="null", Double_t scale2=0.5, TString align="R") {
  // Making -1 a placeholder for function's default value
  if (xpos == -1) xpos = 0.9;
  if (ypos == -1) ypos = 0.9;
  if (scale == -1) scale = 1;
  if (str == "-1") str = "null";
  if (scale2 == -1) scale2 = 0.5;
  if (align == "-1") align = "R";
  // -2 as the first parameter triggers printing comamnd options
  if (xpos == -2) {
    cout
    << endl
    << "  USAGE: BABARSmartLabel(xpos,ypos,scale,\"str\",scale2,\"align\");" << endl
    << "  Prints the official BaBar label on the active ROOT pad" << endl
    << "      xpos    X position of the \"BaBar\" label\'s top right corner, 0 < xpos < 1, defaults to 0.9" << endl
    << "      ypos    Y position of the \"BaBar\" label\'s top right corner, 0 < ypos < 1, defaults to 0.9" << endl
    << "      scale   relative size of the label, defaults to 1" << endl
    << "      str     LaTeX-style text that goes under the BaBar label. Use # instead of \\" << endl
    << "      scale2  relative size of the second line of text, defaults to 0.5" << endl
    << "      align   R or L: str is aligned to the right (default) or left edge of the \"BaBar\" label" << endl
    << "    By default, the second line of text is printed in Helvetica (a sans-serif font)." << endl
    << "  You can use the #font[] command to change the font. Refer to ROOT documentation for" << endl
    << "  more information on use of text in ROOT." << endl;
    cout
    << "    \"Magic\" options: " << endl
    << "      xpos = -2 displays BABARSmartLabel() help" << endl
    << "      -1 can be used as a place holder to use the default value of any of the parameters" << endl
    << "      There are a few predefined values of str that start with a ~ (tilde):" << endl
    << "        \"~1\"       =          \"preliminary\"      " << endl
    << "        \"~2\"       =          \"very preliminary\" " << endl
    << "        \"~2000\"    =          \"year 2000 preliminary\" " << endl
    << "        \"~2001\"    =          \"year 2001 preliminary\" " << endl
    << "        \"~25\"      =          \"25 fb^{-1}\" " << endl
    << "        \"~B->fcp\"  =           a big formula you should have seen before" << endl
    << "    Examples: " << endl
    << "      BABARSmartLabel(); " << endl
    << "      BABARSmartLabel(-1,-1,-1,\"preliminary\"); " << endl
    << "      BABARSmartLabel(0.9,0.8,1.2,\"~25\"); " << endl
    << "      BABARSmartLabel(0.9,0.7,1.2,\"25 fb^{-1} preliminary\",0.6); " << endl
    << "      BABARSmartLabel(0.9,0.5,-1,\"~B->fcp\",0.25); " << endl
    << endl;
    return();
  }
  // A few predefined labels to go to the second line of text
  if (str == "~1") str = "preliminary";
  if (str == "~2") str = "very preliminary";
  if (str == "~2000") str = "year 2000 preliminary";
  if (str == "~2001") str = "year 2000 preliminary";
  if (str == "~25") str = "25 fb^{-1}";
  if (str == "~B->fcp") str = "#font[12]{#Gamma#font[132]{(}B^{#font[12]{0}}_{#font[132]{phys}}#font[132]{(}t#font[132]{)} #rightarrow f_{CP}#font[132]{)} = #left|A_{f_{CP}}#right|^{#font[132]{2}} e^{-#Gamma t} #left[#frac{1+|#lambda_{f_{CP}}|^{#font[132]{2}}}{#font[132]{2}} + #frac{#font[132]{1}-|#lambda_{f_{CP}}|^{#font[132]{2}}}{#font[132]{2}} #font[132]{cos}#font[132]{(}#DeltaMt#font[132]{)} - #font[132]{Im }#lambda_{f_{CP}}#font[132]{sin}#font[132]{(}#DeltaMt#font[132]{)}#right]}";
  
  // Draw the label
  TLatex *babar = new TLatex();
  Double_t cheburashkaFactorX=1, cheburashkaFactorY=1, padSizeX=500, padSizeY=500, xpos2, ypos2, xposL;
  babar->SetNDC(kTRUE);
  babar->SetTextFont(32); // Bold-Italic Times
  babar->SetTextAlign(31); // Right-Bottom
  padSizeX = gPad->GetWw()*gPad->GetWNDC(); // Get pad's dimensions
  padSizeY = gPad->GetWh()*gPad->GetHNDC();
  if (padSizeX>padSizeY) cheburashkaFactorX=padSizeY/padSizeX;
  if (padSizeX<padSizeY) cheburashkaFactorY=padSizeX/padSizeY;
  //xpos2=xpos-0.185*scale*cheburashkaFactorX;
  xpos2=xpos-0.188*scale*cheburashkaFactorX;
  ypos2=ypos-0.0620*scale*cheburashkaFactorY;
  xposL=xpos-0.253*scale*cheburashkaFactorX;
  babar->SetTextSize(0.10*scale); // Beginning to draw "BaBar"
  babar->DrawText(xpos2,ypos2,"G ");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.059*scale*cheburashkaFactorX,ypos2,"12 ");
  babar->SetTextSize(0.10*scale);
  babar->DrawText(xpos2+0.1015*scale*cheburashkaFactorX,ypos2,"J");
  babar->SetTextSize(0.075*scale);
  babar->DrawText(xpos2+0.2275*scale*cheburashkaFactorX,ypos2,"LAB");
  if (str == "null") return();    // Beginning to draw the second line of text
  babar->SetTextFont(42); // Helvetica (medium, upright)
  babar->SetTextSize(0.1*scale2);
  if (align == "L") then {
    babar->SetTextAlign(13); // Left-Top
    babar->DrawLatex(xposL,ypos2-0.02*scale2*cheburashkaFactorY,str);
  }
  else {
    babar->SetTextAlign(33); // Right-Top
    babar->DrawLatex(xpos,ypos2-0.02*scale2*cheburashkaFactorY,str);
  }
  delete babar;
}
