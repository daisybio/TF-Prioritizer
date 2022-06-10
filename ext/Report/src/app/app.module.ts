import {NgModule} from '@angular/core';
import {BrowserModule} from '@angular/platform-browser';

import {AppRoutingModule} from "./app-routing";
import {AppComponent} from './app.component';
import {HeaderComponent} from './templates/accordion/header/header.component';
import {TfGroupAccordionComponent} from './routes/tf-ranking/tf-group-accordion/tf-group-accordion.component';
import {SingleTfAccordionComponent} from './routes/tf-ranking/single-tf-accordion/single-tf-accordion.component';
import {TfRankingComponent} from './routes/tf-ranking/tf-ranking.component';
import {DataHeaderComponent} from './routes/tf-ranking/single-tf-accordion/data-header/data-header.component';
import {DataContentComponent} from './routes/tf-ranking/single-tf-accordion/data-content/data-content.component';
import {ValidationComponent} from './routes/validation/validation.component';
import {DistributionComponent} from './routes/distribution/distribution.component';
import {RegressionComponent} from './routes/regression/regression.component';
import {DocumentationComponent} from './routes/documentation/documentation.component';
import {ParametersComponent} from './routes/parameters/parameters.component';
import {ImageSelectorComponent} from './templates/accordion/image-selector/image-selector.component';
import {HashLocationStrategy, LocationStrategy} from "@angular/common";
import {CrossLinksComponent} from './templates/accordion/cross-links/cross-links.component';
import {ImportantLociComponent} from './routes/important-loci/important-loci.component';
import {TopLog2fcComponent} from './routes/top-log2fc/top-log2fc.component';
import {CoOccurrenceAnalysisComponent} from './routes/co-occurrence-analysis/co-occurrence-analysis.component';

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
