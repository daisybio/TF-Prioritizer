import {NgModule} from '@angular/core';
import {BrowserModule} from '@angular/platform-browser';
import {BrowserAnimationsModule} from "@angular/platform-browser/animations";

import {AppComponent} from './app.component';
import {TfRankingComponent} from './components/tf-ranking/tf-ranking.component';
import {HttpClientModule} from "@angular/common/http";
import {MatChipsModule} from "@angular/material/chips";
import {MatExpansionModule} from "@angular/material/expansion";
import {TfDataComponent} from './components/tf-data/tf-data.component';
import {NgOptimizedImage} from "@angular/common";
import {ImageSelectorComponent} from "./components/image-selector/image-selector.component";
import {DataSelectorComponent} from "./components/data-selector/data-selector.component";
import {MatFormFieldModule} from "@angular/material/form-field";
import {MatSelectModule} from "@angular/material/select";
import { GeneExpressionComponent } from './components/gene-expression/gene-expression.component';

@NgModule({
  declarations: [
    AppComponent,
    TfRankingComponent,
    TfDataComponent,
    ImageSelectorComponent,
    DataSelectorComponent,
    GeneExpressionComponent
  ],
  imports: [
    BrowserModule,
    HttpClientModule,
    MatChipsModule,
    MatExpansionModule,
    BrowserAnimationsModule,
    NgOptimizedImage,
    MatFormFieldModule,
    MatSelectModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule {
}
