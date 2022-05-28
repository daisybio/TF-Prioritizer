import {NgModule} from '@angular/core';
import {BrowserModule} from '@angular/platform-browser';

import {AppRoutingModule} from "./app-routing";
import {AppComponent} from './app.component';
import {HeaderComponent} from './accordion/header/header.component';
import {TfGroupAccordionComponent} from './tf-ranking/tf-group-accordion/tf-group-accordion.component';
import {SingleTfAccordionComponent} from './tf-ranking/single-tf-accordion/single-tf-accordion.component';
import {TfRankingComponent} from './tf-ranking/tf-ranking.component';
import {DataHeaderComponent} from './tf-ranking/single-tf-accordion/data-header/data-header.component';
import {DataContentComponent} from './tf-ranking/single-tf-accordion/data-content/data-content.component';
import {ValidationComponent} from './validation/validation.component';
import {DistributionComponent} from './distribution/distribution.component';
import {RegressionComponent} from './regression/regression.component';
import {DocumentationComponent} from './documentation/documentation.component';
import {ParametersComponent} from './parameters/parameters.component';
import {ImageSelectorComponent} from './accordion/image-selector/image-selector.component';
import {HashLocationStrategy, LocationStrategy} from "@angular/common";
import {CrossLinksComponent} from './accordion/cross-links/cross-links.component';
import {ImportantLociComponent} from './important-loci/important-loci.component';
import {TopLog2fcComponent} from './top-log2fc/top-log2fc.component';
import {CoOccurrenceAnalysisComponent} from './co-occurrence-analysis/co-occurrence-analysis.component';

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    TfGroupAccordionComponent,
    SingleTfAccordionComponent,
    TfRankingComponent,
    DataHeaderComponent,
    DataContentComponent,
    ValidationComponent,
    DistributionComponent,
    RegressionComponent,
    DocumentationComponent,
    ParametersComponent,
    ImageSelectorComponent,
    CrossLinksComponent,
    ImportantLociComponent,
    TopLog2fcComponent,
    CoOccurrenceAnalysisComponent,
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
