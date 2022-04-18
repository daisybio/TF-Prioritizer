import {NgModule} from '@angular/core';
import {BrowserModule} from '@angular/platform-browser';

import {AppRoutingModule} from "./app-routing";
import {AppComponent} from './app.component';
import {HeaderComponent} from './accordion/header/header.component';
import {TfGroupAccordionComponent} from './tf-ranking/tf-group-accordion/tf-group-accordion.component';
import {SingleTfAccordionComponent} from './tf-ranking/single-tf-accordion/single-tf-accordion.component';
import {TfRankingComponent} from './tf-ranking/tf-ranking.component';
import {DataHeaderComponent} from './tf-ranking/single-tf-accordion/data-header/data-header.component';
import {GeneIDComponent} from './tf-ranking/single-tf-accordion/gene-id/gene-id.component';
import {DataContentComponent} from './tf-ranking/single-tf-accordion/data-content/data-content.component';
import {ValidationComponent} from './validation/validation.component';
import {DistributionComponent} from './distribution/distribution.component';
import {RegressionComponent} from './regression/regression.component';
import {DocumentationComponent} from './documentation/documentation.component';
import {ParametersComponent} from './parameters/parameters.component';
import {ImageSelectorComponent} from './accordion/image-selector/image-selector.component';
import {HashLocationStrategy, LocationStrategy} from "@angular/common";

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    TfGroupAccordionComponent,
    SingleTfAccordionComponent,
    TfRankingComponent,
    DataHeaderComponent,
    GeneIDComponent,
    DataContentComponent,
    ValidationComponent,
    DistributionComponent,
    RegressionComponent,
    DocumentationComponent,
    ParametersComponent,
    ImageSelectorComponent,
  ],
  imports: [
    BrowserModule,
    AppRoutingModule
  ],
  providers: [
    {provide: LocationStrategy, useClass: HashLocationStrategy}
  ],
  bootstrap: [AppComponent]
})
export class AppModule {
}
