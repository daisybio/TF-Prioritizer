import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';

import { AppComponent } from './app.component';
import { HeaderComponent } from './accordion/header/header.component';
import { TfGroupAccordionComponent } from './accordion/tf-group-accordion/tf-group-accordion.component';
import { SingleTfAccordionComponent } from './accordion/single-tf-accordion/single-tf-accordion.component';
import { TfGroupsContainerComponent } from './accordion/tf-groups-container/tf-groups-container.component';

@NgModule({
  declarations: [
    AppComponent,
    HeaderComponent,
    TfGroupAccordionComponent,
    SingleTfAccordionComponent,
    TfGroupsContainerComponent,
  ],
  imports: [
    BrowserModule
  ],
  providers: [],
  bootstrap: [AppComponent]
})
export class AppModule { }
